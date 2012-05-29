package lu.jsuper;

import lu.jsuper.Dlu_pdsp_defs.GlobalLU_t;
import lu.jsuper.Dlu_pdsp_defs.pdgstrf_threadarg_t;
import lu.jsuper.Dlu_pdsp_defs.pxgstrf_shared_t;
import lu.jsuper.Dlu_slu_mt_util.Gstat_t;
import lu.jsuper.Dlu_slu_mt_util.superlumt_options_t;
import lu.jsuper.Dlu_slu_mt_util.yes_no_t;
import lu.jsuper.Dlu_supermatrix.SuperMatrix;

import static lu.jsuper.Dlu.PROFILE;
import static lu.jsuper.Dlu.PREDICT_OPT;
import static lu.jsuper.Dlu.DEBUGlevel;
import static lu.jsuper.Dlu.printf;

import static lu.jsuper.Dlu_slu_mt_util.EMPTY;
import static lu.jsuper.Dlu_slu_mt_util.NO_MARKER;
import static lu.jsuper.Dlu_slu_mt_util.TIC;
import static lu.jsuper.Dlu_slu_mt_util.LOCOL;
import static lu.jsuper.Dlu_slu_mt_util.HICOL;
import static lu.jsuper.Dlu_slu_mt_util.BADPAN;
import static lu.jsuper.Dlu_slu_mt_util.BADCOL;
import static lu.jsuper.Dlu_util.ifill;

import static lu.jsuper.Dlu_superlu_timer.SuperLU_timer_;

import static lu.jsuper.Dlu.fflush;
import static lu.jsuper.Dlu.stdout;


public class Dlu_pdgstrf_thread {

	static
	Object[]
	pdgstrf_thread(Object arg[])
	{
	/*
	 * -- SuperLU MT routine (version 2.0) --
	 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
	 * and Xerox Palo Alto Research Center.
	 * September 10, 2007
	 *
	 *
	 * Purpose
	 * =======
	 *
	 * This is the slave process, representing the main scheduling loop to
	 * perform the factorization. Each process executes a copy of the
	 * following code ... (SPMD paradigm)
	 *
	 * Working arrays local to each process
	 * ======================================
	 *   marker[0:3*m-1]: marker[i] == j means node i has been reached when
	 *                                 working on column j.
	 *	Storage: relative to original row subscripts
	 *
	 *	THERE ARE 3 OF THEM:
	 *          marker[0 : m-1]:   used by pdgstrf_factor_snode() and
	 *                                     pdgstrf_panel_dfs();
	 *          marker[m : 2m-1]:  used by pdgstrf_panel_dfs() and
	 *                                     pxgstrf_super_bnd_dfs();
	 *                values in [0 : n-1]  when used by pdgstrf_panel_dfs()
	 *                values in [n : 2n-1] when used by pxgstrf_super_bnd_dfs()
	 *	    marker[2m : 3m-1]: used by pdgstrf_column_dfs() in inner-factor
	 *
	 *   parent[0:n-1]: parent vector used during dfs
	 *      Storage: relative to new row subscripts
	 *
	 *   xplore[0:2m-1]: xplore[i] gives the location of the next (dfs)
	 *	unexplored neighbor of i in lsub[*]; xplore[n+i] gives the
	 *      location of the last unexplored neighbor of i in lsub[*].
	 *
	 *   segrep[0:nseg-1]: contains the list of supernodal representatives
	 *	in topological order of the dfs. A supernode representative is the
	 *	last column of a supernode.
	 *
	 *   repfnz[0:m-1]: for a nonzero segment U[*,j] that ends at a
	 *	supernodal representative r, repfnz[r] is the location of the first
	 *	nonzero in this segment.  It is also used during the dfs:
	 *      repfnz[r]>0 indicates that supernode r has been explored.
	 *	NOTE: There are w of them, each used for one column of a panel.
	 *
	 *   panel_lsub[0:w*m-1]: temporary for the nonzero row indices below
	 *      the panel diagonal. These are filled in during pdgstrf_panel_dfs(),
	 *      and are used later in the inner LU factorization.
	 *	panel_lsub[]/dense[] pair forms the SPA data structure.
	 *	NOTE: There are w of them.
	 *
	 *   dense[0:w*m-1]: sparse accumulator (SPA) for intermediate values;
	 *	NOTE: there are w of them.
	 *
	 *   tempv[0:m-1]: real temporary used for dense numeric kernels;
	 *
	 *
	 * Scheduling algorithm (For each process ...)
	 * ====================
	 *     Shared task Q <-- { relaxed s-nodes (CANGO) };
	 *
	 *     WHILE (not finished)
	 *
	 *         panel = Scheduler(Q); (see pxgstrf_scheduler.c for policy)
	 *
	 *         IF (panel == RELAXED_SNODE)
	 *             factor_relax_snode(panel);
	 *         ELSE
	 *             * pdgstrf_panel_dfs()
	 *                 - skip all BUSY s-nodes (or panels)
	 *
	 *             * dpanel_bmod()
	 *                 - updates from DONE s-nodes
	 *                 - wait for BUSY s-nodes to become DONE
	 *
	 *             * inner-factor()
	 *                 - identical as it is in the sequential algorithm,
	 *                   except that pruning() will interact with the
	 *                   pdgstrf_panel_dfs() of other panels.
	 *         ENDIF
	 *
	 *     END WHILE;
	 *
	 */

	    pdgstrf_threadarg_t thr_arg = arg;
	    int         pnum = thr_arg.pnum;

	    /* Unpack the options argument */
	    superlumt_options_t superlumt_options = thr_arg.superlumt_options;
	    pxgstrf_shared_t  pxgstrf_shared = thr_arg.pxgstrf_shared;
	    int         panel_size = superlumt_options.panel_size;
	    double     diag_pivot_thresh = superlumt_options.diag_pivot_thresh;
	    /* may be modified */
	    yes_no_t    usepr[]     = superlumt_options.usepr;
	    int         etree[]     = superlumt_options.etree;
	    int         super_bnd[] = superlumt_options.part_super_h;
	    int         perm_r[]    = superlumt_options.perm_r;
	    int         inv_perm_c[]= pxgstrf_shared.inv_perm_c;
	    int         inv_perm_r[]= pxgstrf_shared.inv_perm_r;
	    int	        xprune[]    = pxgstrf_shared.xprune;
	    int	        ispruned[]  = pxgstrf_shared.ispruned;
	    SuperMatrix A          = pxgstrf_shared.A;
	    GlobalLU_t  Glu        = pxgstrf_shared.Glu;
	    Gstat_t 	Gstat      = pxgstrf_shared.Gstat;
	    int         info[]     = thr_arg.info;

	    /* Local working arrays */
	    int       iwork[];
	    double    dwork[];
	    int	      segrep[], repfnz[], parent[], xplore[];
	    int	      panel_lsub[]; /* dense[]/panel_lsub[] pair forms a w-wide SPA */
	    int	      marker[], marker1[], marker2[];
	    int       lbusy[]; /* "Local busy" array, indicates which descendants
				 were busy when this panel's computation began.
				 Those columns (s-nodes) are treated specially
				 during pdgstrf_panel_dfs() and dpanel_bmod(). */

	    int       spa_marker[]; /* size n-by-w */
	    int       w_lsub_end[]; /* record the end of each column in panel_lsub */
	    double    dense[], tempv[];
	    int       lsub[], xlsub[], xlsub_end[];

	    /* Local scalars */
	    int m, n, k, jj, jcolm1, itemp, singular;
	    int       pivrow;   /* pivotal row number in the original matrix A */
	    int       nseg1;	/* no of segments in U-column above panel row jcol */
	    int       nseg;	/* no of segments in each U-column */
	    int       w, bcol, jcol;

	if (PROFILE) {
	    double utime[] = Gstat.utime;
	    double t1, t2, t, stime;
	    float flopcnt;
	}

	if (PREDICT_OPT) {
	    flops_t  ops[] = Gstat.ops;
	    float pdiv;
	}

	if ( DEBUGlevel>=1 ) {
	    printf("(%d) thr_arg. pnum %d, info %d\n", pnum, thr_arg.pnum, thr_arg.info);
	}

	    singular   = 0;
	    m          = A.nrow;
	    n          = A.ncol;
	    lsub       = Glu.lsub;
	    xlsub      = Glu.xlsub;
	    xlsub_end  = Glu.xlsub_end;

	    /* Allocate and initialize the per-process working storage. */
	    if ( (info[0] = pdgstrf_WorkInit(m, panel_size, &iwork, &dwork)) ) {
		*info += pdgstrf_memory_use(Glu.nzlmax, Glu.nzumax, Glu.nzlumax);
		return 0;
	    }
	    pxgstrf_SetIWork(m, panel_size, iwork, &segrep, &parent, &xplore,
		     &repfnz, &panel_lsub, &marker, &lbusy);
	    pdgstrf_SetRWork(m, panel_size, dwork, &dense, &tempv);

	    /* New data structures to facilitate parallel algorithm */
	    spa_marker = intMalloc(m * panel_size);
	    w_lsub_end = intMalloc(panel_size);
	    ifill (spa_marker, m * panel_size, EMPTY);
	    ifill (marker, m * NO_MARKER, EMPTY);
	    ifill (lbusy, m, EMPTY);
	    jcol = EMPTY;
	    marker1 = marker + m;
	    marker2 = marker + 2*m;

	if (PROFILE) {
	    stime = SuperLU_timer_();
	}

	    /* -------------------------
	       Main loop: repeatedly ...
	       ------------------------- */
	    while ( pxgstrf_shared.tasks_remain > 0 ) {

	if (PROFILE) {
		TIC(t);
	}
		/* Get a panel from the scheduler. */
		pxgstrf_scheduler(pnum, n, etree, &jcol, &bcol, pxgstrf_shared);

	if ( DEBUGlevel>=1 ) {
	    if ( jcol>=LOCOL && jcol<=HICOL ) {
		printf("(%d) Scheduler(): jcol %d, bcol %d, tasks_remain %d\n",
		       pnum, jcol, bcol, pxgstrf_shared.tasks_remain);
		fflush(stdout);
	    }
	}

	if (PROFILE) {
		TOC(t2, t);
		Gstat.procstat[pnum].skedtime += t2;
	}

		if ( jcol != EMPTY ) {
		    w = pxgstrf_shared.pan_status[jcol].size;

	if ( DEBUGlevel>=3 ) {
		    printf("P%2d got panel %5d-%5d\ttime %.4f\tpanels_left %d\n",
			   pnum, jcol, jcol+w-1, SuperLU_timer_(),
			   pxgstrf_shared.tasks_remain);
		    fflush(stdout);
	}
		    /* Nondomain panels */
	if (PROFILE) {
		    flopcnt = Gstat.procstat[pnum].fcops;
		    Gstat.panstat[jcol].pnum = pnum;
		    TIC(t1);
		    Gstat.panstat[jcol].starttime = t1;
	}
		    if ( pxgstrf_shared.pan_status[jcol].type == RELAXED_SNODE ) {

	if (PREDICT_OPT) {
			pdiv = Gstat.procstat[pnum].fcops;
	}
			/* A relaxed supernode at the bottom of the etree */
			pdgstrf_factor_snode
			    (pnum, jcol, A, diag_pivot_thresh, usepr,
			     perm_r, inv_perm_r, inv_perm_c, xprune, marker,
			     panel_lsub, dense, tempv, pxgstrf_shared, info);
			if ( info[0] ) {
			    if ( info[0] > n ) return 0;
			    else if ( singular == 0 || info[0] < singular )
			        singular = info[0];
	if ( DEBUGlevel>=1 ) {
	    printf("(%d) After pdgstrf_factor_snode(): singular=%d\n", pnum, singular);
	}
			}

			/* Release the whole relaxed supernode */
			for (jj = jcol; jj < jcol + w; ++jj)
			    pxgstrf_shared.spin_locks[jj] = 0;
	if (PREDICT_OPT) {
			pdiv = Gstat.procstat[pnum].fcops - pdiv;
			cp_panel[jcol].pdiv = pdiv;
	}
		    } else { /* Regular panel */
	if (PROFILE) {
			TIC(t);
	}
			pxgstrf_mark_busy_descends(pnum, jcol, etree, pxgstrf_shared,
						   &bcol, lbusy);

			/* Symbolic factor on a panel of columns */
			pdgstrf_panel_dfs
			    (pnum, m, w, jcol, A, perm_r, xprune,ispruned,lbusy,
			     &nseg1, panel_lsub, w_lsub_end, segrep, repfnz,
			     marker, spa_marker, parent, xplore, dense, Glu);
	if ( DEBUGlevel>=2 ) {
	  if ( jcol==BADPAN )
	    printf("(%d) After pdgstrf_panel_dfs(): nseg1 %d, w_lsub_end %d\n",
		   pnum, nseg1, w_lsub_end[0]);
	}
	if (PROFILE) {
			TOC(t2, t);
			utime[DFS] += t2;
	}
			/* Numeric sup-panel updates in topological order.
			 * On return, the update values are temporarily stored in
			 * the n-by-w SPA dense[m,w].
			 */
			pdgstrf_panel_bmod
			    (pnum, m, w, jcol, bcol, inv_perm_r, etree,
			     &nseg1, segrep, repfnz, panel_lsub, w_lsub_end,
			     spa_marker, dense, tempv, pxgstrf_shared);

			/*
			 * All "busy" descendants are "done" now --
			 * Find the set of row subscripts in the preceeding column
			 * "jcol-1" of the current panel. Column "jcol-1" is
			 * usually taken by a process other than myself.
			 * This row-subscripts information will be used by myself
			 * during column dfs to detect whether "jcol" belongs
			 * to the same supernode as "jcol-1".
			 *
			 * ACCORDING TO PROFILE, THE AMOUNT OF TIME SPENT HERE
			 * IS NEGLIGIBLE.
			 */
			jcolm1 = jcol - 1;
			itemp = xlsub_end[jcolm1];
			for (k = xlsub[jcolm1]; k < itemp; ++k)
			    marker2[lsub[k]] = jcolm1;
	if (PREDICT_OPT) {
			pdiv = Gstat.procstat[pnum].fcops;
	}
			/* Inner-factorization, using sup-col algorithm */
			for ( jj = jcol; jj < jcol + w; jj++) {
			    k = (jj - jcol) * m; /* index into w-wide arrays */
			    nseg = nseg1; /* begin after all the panel segments */
	if (PROFILE) {
			    TIC(t);
	}
			    /* Allocate storage for the current H-supernode. */
			    if ( Glu.dynamic_snode_bound != 0 && super_bnd[jj] != 0 ) {
			        /* jj starts a supernode in H */
				pxgstrf_super_bnd_dfs
				    (pnum, m, n, jj, super_bnd[jj], A, perm_r,
				     inv_perm_r, xprune, ispruned, marker1, parent,
				     xplore, pxgstrf_shared);
			    }

			    if ( (*info = pdgstrf_column_dfs
				            (pnum, m, jj, jcol, perm_r, ispruned,
					     &panel_lsub[k],w_lsub_end[jj-jcol],
					     super_bnd, &nseg, segrep,
					     &repfnz[k], xprune, marker2,
					     parent, xplore, pxgstrf_shared)) )
				return 0;
	if (PROFILE) {
			    TOC(t2, t);
			    utime[DFS] += t2;
	}
			    /* On return, the L supernode is gathered into the
			       global storage. */
			    if ( (info[0] = pdgstrf_column_bmod
				          (pnum, jj, jcol, (nseg - nseg1),
					   &segrep[nseg1], &repfnz[k],
					   &dense[k], tempv, pxgstrf_shared, Gstat)) )
				return 0;

			    if ( (info[0] = pdgstrf_pivotL
				            (pnum, jj, diag_pivot_thresh, usepr,
					     perm_r, inv_perm_r, inv_perm_c,
					     &pivrow, Glu, Gstat)) )
				if ( singular == 0 || info[0] < singular ) {
				    singular = info[0];
	if ( DEBUGlevel>=1 ) {
	    printf("(%d) After pdgstrf_pivotL(): singular=%d\n", pnum, singular);
	}
				}

	                    /* release column "jj", so that the other processes
	                       waiting for this column can proceed */
			    pxgstrf_shared.spin_locks[jj] = 0;

			    /* copy the U-segments to ucol[*] */
			    if ( (info[0] = pdgstrf_copy_to_ucol
				            (pnum,jj,nseg,segrep,&repfnz[k],
					     perm_r, &dense[k], pxgstrf_shared)) )
			      return 0;

			    /* Prune columns [0:jj-1] using column jj */
			    pxgstrf_pruneL(jj, perm_r, pivrow, nseg, segrep,
					   &repfnz[k], xprune, ispruned, Glu);

			    /* Reset repfnz[] for this column */
			    pxgstrf_resetrep_col (nseg, segrep, &repfnz[k]);

	if ( DEBUGlevel>=2 ) {
	/*  if (jj >= LOCOL && jj <= HICOL) {*/
	  if ( jj==BADCOL ) {
	    dprint_lu_col(pnum, "panel:", jcol, jj, w, pivrow, xprune, Glu);
	    dcheck_zero_vec(pnum, "after pdgstrf_copy_to_ucol() dense_col[]", n, &dense[k]);
	  }
	}
			} /* for jj ... */

	if (PREDICT_OPT) {
			pdiv = Gstat.procstat[pnum].fcops - pdiv;
			cp_panel[jcol].pdiv = pdiv;
	}

		    } /* else regular panel ... */

		    STATE( jcol ) = DONE; /* Release panel jcol. */

	if (PROFILE) {
		    TOC(Gstat.panstat[jcol].fctime, t1);
		    Gstat.panstat[jcol].flopcnt += Gstat.procstat[pnum].fcops - flopcnt;
		    /*if ( Glu.tasks_remain < P ) {
			flops_last_P_panels += Gstat.panstat[jcol].flopcnt;
			printf("Panel %d, flops %e\n", jcol, Gstat.panstat[jcol].flopcnt);
			fflush(stdout);
		    } */
	}

		}
	if (PROFILE) {
		else { /* No panel from the task queue - wait and try again */
		    Gstat.procstat[pnum].skedwaits++;
		}
	}

	    } /* while there are more panels */

	    info[0] = singular;

	    /* Free work space and compress storage */
	    pdgstrf_WorkFree(iwork, dwork, Glu);

	if (PROFILE) {
	    Gstat.procstat[pnum].fctime = SuperLU_timer_() - stime;
	}

	    return 0;
	}


}

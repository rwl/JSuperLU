package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_pdsp_defs.pdgstrf_threadarg_t;
import gov.lbl.superlu.Dlu_pdsp_defs.pxgstrf_shared_t;
import gov.lbl.superlu.Dlu_slu_mt_util.Gstat_t;
import gov.lbl.superlu.Dlu_slu_mt_util.superlumt_options_t;
import gov.lbl.superlu.Dlu_slu_mt_util.yes_no_t;
import gov.lbl.superlu.Dlu_supermatrix.SuperMatrix;

import static gov.lbl.superlu.Dlu.PROFILE;
import static gov.lbl.superlu.Dlu.PREDICT_OPT;
import static gov.lbl.superlu.Dlu.DEBUGlevel;
import static gov.lbl.superlu.Dlu.printf;

import static gov.lbl.superlu.Dlu_slu_mt_util.EMPTY;
import static gov.lbl.superlu.Dlu_slu_mt_util.NO_MARKER;
import static gov.lbl.superlu.Dlu_slu_mt_util.TIC;
import static gov.lbl.superlu.Dlu_slu_mt_util.TOC;
import static gov.lbl.superlu.Dlu_slu_mt_util.LOCOL;
import static gov.lbl.superlu.Dlu_slu_mt_util.HICOL;
import static gov.lbl.superlu.Dlu_slu_mt_util.BADPAN;
import static gov.lbl.superlu.Dlu_slu_mt_util.BADCOL;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.DFS;

import static gov.lbl.superlu.Dlu_util.ifill;
import static gov.lbl.superlu.Dlu_util.pxgstrf_resetrep_col;

import static gov.lbl.superlu.Dlu_pdutil.dprint_lu_col;
import static gov.lbl.superlu.Dlu_pdutil.dcheck_zero_vec;

import static gov.lbl.superlu.Dlu_superlu_timer.SuperLU_timer_;

import static gov.lbl.superlu.Dlu.fflush;
import static gov.lbl.superlu.Dlu.stdout;

import static gov.lbl.superlu.Dlu_pdmemory.pdgstrf_WorkInit;
import static gov.lbl.superlu.Dlu_pdmemory.pdgstrf_memory_use;
import static gov.lbl.superlu.Dlu_pdmemory.pdgstrf_SetRWork;
import static gov.lbl.superlu.Dlu_pdmemory.pdgstrf_WorkFree;

import static gov.lbl.superlu.Dlu_pmemory.pxgstrf_SetIWork;
import static gov.lbl.superlu.Dlu_pmemory.intMalloc;

import static gov.lbl.superlu.Dlu_pxgstrf_scheduler.pxgstrf_scheduler;
import static gov.lbl.superlu.Dlu_pxgstrf_synch.panel_t.RELAXED_SNODE;
import static gov.lbl.superlu.Dlu_pxgstrf_synch.pipe_state_t.DONE;
import static gov.lbl.superlu.Dlu_pdgstrf_factor_snode.pdgstrf_factor_snode;
import static gov.lbl.superlu.Dlu_pxgstrf_mark_busy_descends.pxgstrf_mark_busy_descends;
import static gov.lbl.superlu.Dlu_pdgstrf_panel_dfs.pdgstrf_panel_dfs;
import static gov.lbl.superlu.Dlu_pdgstrf_panel_bmod.pdgstrf_panel_bmod;
import static gov.lbl.superlu.Dlu_pxgstrf_super_bnd_dfs.pxgstrf_super_bnd_dfs;
import static gov.lbl.superlu.Dlu_pdgstrf_column_dfs.pdgstrf_column_dfs;
import static gov.lbl.superlu.Dlu_pdgstrf_column_bmod.pdgstrf_column_bmod;
import static gov.lbl.superlu.Dlu_pdgstrf_pivotL.pdgstrf_pivotL;
import static gov.lbl.superlu.Dlu_pdgstrf_copy_to_ucol.pdgstrf_copy_to_ucol;
import static gov.lbl.superlu.Dlu_pxgstrf_pruneL.pxgstrf_pruneL;


public class Dlu_pdgstrf_thread {

	static
	Object[]
	pdgstrf_thread(pdgstrf_threadarg_t arg)
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
	    int       iwork[][] = new int[1][];
	    double    dwork[][] = new double[1][];
	    int[][]	      segrep, repfnz, parent, xplore;
	    segrep = new int[1][];
	    repfnz = new int[1][];
	    parent = new int[1][];
	    xplore = new int[1][];
	    int	      panel_lsub[][]; /* dense[]/panel_lsub[] pair forms a w-wide SPA */
	    panel_lsub = new int[1][];
	    int	      marker[][], marker1[], marker2[];
	    int marker1_offset, marker2_offset;
	    marker = new int[1][];
	    int       lbusy[][]; /* "Local busy" array, indicates which descendants
				 were busy when this panel's computation began.
				 Those columns (s-nodes) are treated specially
				 during pdgstrf_panel_dfs() and dpanel_bmod(). */
	    lbusy = new int[1][];

	    int       spa_marker[]; /* size n-by-w */
	    int       w_lsub_end[]; /* record the end of each column in panel_lsub */
	    double    dense[][] = new double[1][], tempv[][] = new double[1][];
	    int       lsub[], xlsub[], xlsub_end[];

	    /* Local scalars */
	    int m, n, k, jj, jcolm1, itemp, singular;
	    int       pivrow[];   /* pivotal row number in the original matrix A */
	    pivrow = new int[1];
	    int       nseg1[];	/* no of segments in U-column above panel row jcol */
	    nseg1 = new int[1];
	    int       nseg[];	/* no of segments in each U-column */
	    nseg = new int[1];
	    int       w, bcol[], jcol[];
	    bcol = new int[1];
	    jcol = new int[1];

	    double utime[] = Gstat.utime;
	    double t1[], t2[], t[], stime = 0;
	    t = new double[1];
	    t1 = new double[1];
	    t2 = new double[1];
	    float flopcnt = 0;

	    float  ops[] = Gstat.ops;
	    float pdiv = 0;

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
	    if ( (info[0] = pdgstrf_WorkInit(m, panel_size, iwork, dwork)) != 0 ) {
		info[0] += pdgstrf_memory_use(Glu.nzlmax, Glu.nzumax, Glu.nzlumax);
		return null/*0*/;
	    }
	    pxgstrf_SetIWork(m, panel_size, /*iwork, */segrep, parent, xplore,
		     repfnz, panel_lsub, marker, lbusy);
	    pdgstrf_SetRWork(m, panel_size, dwork[0], dense, tempv);

	    /* New data structures to facilitate parallel algorithm */
	    spa_marker = intMalloc(m * panel_size);
	    w_lsub_end = intMalloc(panel_size);
	    ifill (spa_marker, m * panel_size, EMPTY);
	    ifill (marker[0], m * NO_MARKER, EMPTY);
	    ifill (lbusy[0], m, EMPTY);
	    jcol[0] = EMPTY;
	    marker1 = marker[0];
	    marker1_offset = m;
	    marker2 = marker[0];
	    marker2_offset = 2*m;

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
		pxgstrf_scheduler(pnum, n, etree, jcol, bcol, pxgstrf_shared);

	if ( DEBUGlevel>=1 ) {
	    if ( jcol[0]>=LOCOL && jcol[0]<=HICOL ) {
		printf("(%d) Scheduler(): jcol %d, bcol %d, tasks_remain %d\n",
		       pnum, jcol, bcol, pxgstrf_shared.tasks_remain);
		fflush(stdout);
	    }
	}

	if (PROFILE) {
		TOC(t2, t[0]);
		Gstat.procstat[pnum].skedtime += t2[0];
	}

		if ( jcol[0] != EMPTY ) {
		    w = pxgstrf_shared.pan_status[jcol[0]].size;

	if ( DEBUGlevel>=3 ) {
		    printf("P%2d got panel %5d-%5d\ttime %.4f\tpanels_left %d\n",
			   pnum, jcol[0], jcol[0]+w-1, SuperLU_timer_(),
			   pxgstrf_shared.tasks_remain);
		    fflush(stdout);
	}
		    /* Nondomain panels */
	if (PROFILE) {
		    flopcnt = Gstat.procstat[pnum].fcops;
		    Gstat.panstat[jcol[0]].pnum = pnum;
		    TIC(t1);
		    Gstat.panstat[jcol[0]].starttime = t1[0];
	}
		    if ( pxgstrf_shared.pan_status[jcol[0]].type == RELAXED_SNODE ) {

	if (PREDICT_OPT) {
			pdiv = Gstat.procstat[pnum].fcops;
	}
			/* A relaxed supernode at the bottom of the etree */
			pdgstrf_factor_snode
			    (pnum, jcol[0], A, diag_pivot_thresh, usepr,
			     perm_r, inv_perm_r, inv_perm_c, xprune, marker[0],
			     panel_lsub[0], dense[0], tempv[0], pxgstrf_shared, info);
			if ( info[0] != 0 ) {
			    if ( info[0] > n ) return null/*0*/;
			    else if ( singular == 0 || info[0] < singular )
			        singular = info[0];
	if ( DEBUGlevel>=1 ) {
	    printf("(%d) After pdgstrf_factor_snode(): singular=%d\n", pnum, singular);
	}
			}

			/* Release the whole relaxed supernode */
			for (jj = jcol[0]; jj < jcol[0] + w; ++jj)
			    pxgstrf_shared.spin_locks[jj] = 0;
	if (PREDICT_OPT) {
			pdiv = Gstat.procstat[pnum].fcops - pdiv;
			Gstat.cp_panel[jcol[0]].pdiv = pdiv;
	}
		    } else { /* Regular panel */
	if (PROFILE) {
			TIC(t);
	}
			pxgstrf_mark_busy_descends(pnum, jcol[0], etree, pxgstrf_shared,
						   bcol, lbusy[0]);

			/* Symbolic factor on a panel of columns */
			pdgstrf_panel_dfs
			    (pnum, m, w, jcol[0], A, perm_r, xprune,ispruned,lbusy[0],
			     nseg1, panel_lsub[0], w_lsub_end, segrep[0], repfnz[0],
			     marker[0], spa_marker, parent[0], xplore[0], dense[0], Glu);
	if ( DEBUGlevel>=2 ) {
	  if ( jcol[0]==BADPAN )
	    printf("(%d) After pdgstrf_panel_dfs(): nseg1 %d, w_lsub_end %d\n",
		   pnum, nseg1, w_lsub_end[0]);
	}
	if (PROFILE) {
			TOC(t2, t[0]);
			utime[DFS.ordinal()] += t2[0];
	}
			/* Numeric sup-panel updates in topological order.
			 * On return, the update values are temporarily stored in
			 * the n-by-w SPA dense[m,w].
			 */
			pdgstrf_panel_bmod
			    (pnum, m, w, jcol[0], bcol[0], inv_perm_r, etree,
			     nseg1, segrep[0], repfnz[0], panel_lsub[0], w_lsub_end,
			     spa_marker, dense[0], tempv[0], pxgstrf_shared);

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
			jcolm1 = jcol[0] - 1;
			itemp = xlsub_end[jcolm1];
			for (k = xlsub[jcolm1]; k < itemp; ++k)
			    marker2[marker2_offset+lsub[k]] = jcolm1;
	if (PREDICT_OPT) {
			pdiv = Gstat.procstat[pnum].fcops;
	}
			/* Inner-factorization, using sup-col algorithm */
			for ( jj = jcol[0]; jj < jcol[0] + w; jj++) {
			    k = (jj - jcol[0]) * m; /* index into w-wide arrays */
			    nseg[0] = nseg1[0]; /* begin after all the panel segments */
	if (PROFILE) {
			    TIC(t);
	}
			    /* Allocate storage for the current H-supernode. */
			    if ( Glu.dynamic_snode_bound != 0 && super_bnd[jj] != 0 ) {
			        /* jj starts a supernode in H */
				pxgstrf_super_bnd_dfs
				    (pnum, m, n, jj, super_bnd[jj], A, perm_r,
				     inv_perm_r, xprune, ispruned, marker1, marker1_offset, parent[0],
				     xplore[0], pxgstrf_shared);
			    }

			    if ( (info[0] = pdgstrf_column_dfs
				            (pnum, m, jj, jcol[0], perm_r, ispruned,
					     panel_lsub[k],w_lsub_end[jj-jcol[0]],
					     super_bnd, nseg, segrep[0],
					     repfnz[k], xprune, marker2, marker2_offset,
					     parent[0], xplore[0], pxgstrf_shared)) != 0 )
				return null/*0*/;
	if (PROFILE) {
			    TOC(t2, t[0]);
			    utime[DFS.ordinal()] += t2[0];
	}
			    /* On return, the L supernode is gathered into the
			       global storage. */
			    if ( (info[0] = pdgstrf_column_bmod
				          (pnum, jj, jcol[0], (nseg[0] - nseg1[0]),
					   segrep[nseg1[0]], repfnz[k],
					   dense[k], tempv[0], pxgstrf_shared, Gstat)) != 0 )
				return null/*0*/;

			    if ( (info[0] = pdgstrf_pivotL
				            (pnum, jj, diag_pivot_thresh, usepr,
					     perm_r, inv_perm_r, inv_perm_c,
					     pivrow, Glu, Gstat)) != 0 )
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
				            (pnum,jj,nseg[0],segrep[0],repfnz[k],
					     perm_r, dense[k], pxgstrf_shared)) != 0 )
			      return null/*0*/;

			    /* Prune columns [0:jj-1] using column jj */
			    pxgstrf_pruneL(jj, perm_r, pivrow[0], nseg[0], segrep[0],
					   repfnz[k], xprune, ispruned, Glu);

			    /* Reset repfnz[] for this column */
			    pxgstrf_resetrep_col (nseg[0], segrep[0], repfnz[k]);

	if ( DEBUGlevel>=2 ) {
	/*  if (jj >= LOCOL && jj <= HICOL) {*/
	  if ( jj==BADCOL ) {
	    dprint_lu_col(pnum, "panel:", jcol[0], jj, w, pivrow[0], xprune, Glu);
	    dcheck_zero_vec(pnum, "after pdgstrf_copy_to_ucol() dense_col[]", n, dense[k]);
	  }
	}
			} /* for jj ... */

	if (PREDICT_OPT) {
			pdiv = Gstat.procstat[pnum].fcops - pdiv;
			Gstat.cp_panel[jcol[0]].pdiv = pdiv;
	}

		    } /* else regular panel ... */

		    //STATE( jcol[0] ) = DONE; /* Release panel jcol. */
		    pxgstrf_shared.pan_status[jcol[0]].state = DONE;

	if (PROFILE) {
			double[] tx = new double[1];
		    TOC(tx, t1[0]);
		    Gstat.panstat[jcol[0]].fctime = tx[0];
		    Gstat.panstat[jcol[0]].flopcnt += Gstat.procstat[pnum].fcops - flopcnt;
		    /*if ( Glu.tasks_remain < P ) {
			flops_last_P_panels += Gstat.panstat[jcol].flopcnt;
			printf("Panel %d, flops %e\n", jcol, Gstat.panstat[jcol].flopcnt);
			fflush(stdout);
		    } */
	}

		} else {
	if (PROFILE) {
		/* No panel from the task queue - wait and try again */
		Gstat.procstat[pnum].skedwaits++;
	}
		}

	    } /* while there are more panels */

	    info[0] = singular;

	    /* Free work space and compress storage */
	    pdgstrf_WorkFree(iwork[0], dwork[0], Glu);

	if (PROFILE) {
	    Gstat.procstat[pnum].fctime = SuperLU_timer_() - stime;
	}

	    return null/*0*/;
	}

}

package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_pdsp_defs.pxgstrf_shared_t;
import gov.lbl.superlu.Dlu_slu_mt_util.Gstat_t;

import static gov.lbl.superlu.Dlu_sp_ienv.sp_ienv;

import static gov.lbl.superlu.Dlu_slu_mt_util.yes_no_t.YES;
import static gov.lbl.superlu.Dlu_slu_mt_util.yes_no_t.NO;
import static gov.lbl.superlu.Dlu_slu_mt_util.EMPTY;
import static gov.lbl.superlu.Dlu_slu_mt_util.BADCOL;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPER_REP;
import static gov.lbl.superlu.Dlu_slu_mt_util.SINGLETON;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPER_FSUPC;
import static gov.lbl.superlu.Dlu_slu_mt_util.MemType.LSUB;

import static gov.lbl.superlu.Dlu.DEBUGlevel;
import static gov.lbl.superlu.Dlu.printf;
import static gov.lbl.superlu.Dlu.PROFILE;

import static gov.lbl.superlu.Dlu_pxgstrf_synch.NewNsuper;

import static gov.lbl.superlu.Dlu_pmemory.Glu_alloc;
import static gov.lbl.superlu.Dlu_util.PrintInt10;


public class Dlu_pdgstrf_column_dfs {

    static  int  first = 1, maxsuper;

    static
	int
	pdgstrf_column_dfs(
			   final int  pnum,    /* process number */
			   final int  m,       /* number of rows in the matrix */
			   final int  jcol,    /* current column in the panel */
			   final int  fstcol,  /* first column in the panel */
			   int perm_r[],   /* row pivotings that are done so far */
			   int ispruned[], /* in */
			   int col_lsub[], /* the RHS vector to start the dfs */
			   int lsub_end,  /* size of col_lsub[] */
			   int super_bnd[],/* supernode partition by upper bound */
			   int nseg[],     /* modified - with new segments appended */
			   int segrep[],   /* modified - with new segments appended */
			   int repfnz[],   /* modified */
			   int xprune[],   /* modified */
			   int marker2[],  /* modified */
			   int parent[],   /* working array */
			   int xplore[],   /* working array */
			   pxgstrf_shared_t pxgstrf_shared /* modified */
			   )
	{
	/*
	 * -- SuperLU MT routine (version 2.0) --
	 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
	 * and Xerox Palo Alto Research Center.
	 * September 10, 2007
	 *
	 * Purpose
	 * =======
	 *   pdgstrf_column_dfs() performs a symbolic factorization on column jcol,
	 *   and detects whether column jcol belongs in the same supernode as jcol-1.
	 *
	 * Local parameters
	 * ================
	 *   A supernode representative is the last column of a supernode.
	 *   The nonzeros in U[*,j] are segments that end at supernodal
	 *   representatives. The routine returns a list of such supernodal
	 *   representatives in topological order of the dfs that generates them.
	 *   The location of the first nonzero in each such supernodal segment
	 *   (supernodal entry location) is also returned.
	 *
	 *   nseg: no of segments in current U[*,j]
	 *   samesuper: samesuper=NO if column j does not belong in the same
	 *	        supernode as j-1. Otherwise, samesuper=YES.
	 *
	 *   marker2: A-row -. A-row/col (0/1)
	 *   repfnz: SuperA-col -. PA-row
	 *   parent: SuperA-col -. SuperA-col
	 *   xplore: SuperA-col -. index to L-structure
	 *
	 * Return value
	 * ============
	 *     0  success;
	 *   > 0  number of bytes allocated when run out of space.
	 *
	 */
	    GlobalLU_t Glu = pxgstrf_shared.Glu; /* modified */
	    Gstat_t Gstat = pxgstrf_shared.Gstat; /* modified */
	    int jcolm1, jcolm1size, nextl, ifrom;
	    int k, krep, krow, kperm, samesuper, nsuper = 0;
	    int no_lsub;
	    int	    fsupc = 0;		/* first column in a supernode */
	    int     myfnz;		/* first nonz column in a U-segment */
	    int	    chperm, chmark, chrep, kchild;
	    int     xdfs, maxdfs, kpar;
	    int     ito[] = new int[1];	        /* Used to compress row subscripts */
	    int     mem_error;
	    int[]     xsup, xsup_end, supno, lsub, xlsub, xlsub_end;

	    if ( first != 0 ) {
		maxsuper = sp_ienv(3);
		first = 0;
	    }

	    /* Initialize pointers */
	    xsup      = Glu.xsup;
	    xsup_end  = Glu.xsup_end;
	    supno     = Glu.supno;
	    lsub      = Glu.lsub;
	    xlsub     = Glu.xlsub;
	    xlsub_end = Glu.xlsub_end;
	    jcolm1    = jcol - 1;
	    nextl     = lsub_end;
	    no_lsub   = 0;
	    samesuper = YES.ordinal();

	    /* Test whether the row structure of column jcol is contained
	       in that of column jcol-1. */
	    for (k = 0; k < lsub_end; ++k) {
		krow = col_lsub[k];
		if ( perm_r[krow] == EMPTY ) { /* krow is in L */
		    ++no_lsub;
		    if (marker2[krow] != jcolm1)
		        samesuper = NO.ordinal(); /* row subset test */
		    marker2[krow] = jcol;
		}
	    }

	if ( DEBUGlevel>=2 ) {
	  if (jcol == BADCOL)
	    printf("(%d) pdgstrf_column_dfs[1] %d, fstcol %d, lsub_end %d, no_lsub %d, samesuper? %d\n",
		   pnum, jcol, fstcol, lsub_end, no_lsub, samesuper);
	}

	    /*
	     * For each nonzero in A[fstcol:n,jcol] perform DFS ...
	     */
	    for (k = 0; k < lsub_end; ++k) {
		krow = col_lsub[k];

		/* if krow was visited before, go to the next nonzero */
		if ( marker2[krow] == jcol ) continue;
		marker2[krow] = jcol;
		kperm = perm_r[krow];
	if ( DEBUGlevel>=3 ) {
	  if (jcol == BADCOL)
	    printf("(%d) pdgstrf_column_dfs[inner]: perm_r[krow=%d] %d\n", pnum, krow, kperm);
	}

		/* Ignore the nonzeros in U corresponding to the busy columns
		   during the panel DFS. */
		/*if ( lbusy[kperm] != fstcol ) {  xiaoye? */
		if ( kperm >= fstcol ) {
		    /*
		     * krow is in U: if its supernode representative krep
		     * has been explored, update repfnz[*].
		     */
		    krep = SUPER_REP(xsup_end, supno[kperm]);
		    myfnz = repfnz[krep];

	if ( DEBUGlevel>=3 ) {
	  if (jcol == BADCOL)
	    printf("(%d) pdgstrf_column_dfs[inner-U]: krep %d, myfnz %d, kperm %d\n",
		   pnum, krep, myfnz, kperm);
	}
		    if ( myfnz != EMPTY ) {	/* Visited before */
			if ( myfnz > kperm ) repfnz[krep] = kperm;
			/* continue; */
		    } else {
			/* Otherwise, perform dfs starting at krep */
			parent[krep] = EMPTY;
			repfnz[krep] = kperm;
			if ( ispruned[krep] != 0 ) {
			    if ( SINGLETON( xsup_end, xsup_end, supno[krep] ) )
				xdfs = xlsub_end[krep];
			    else xdfs = xlsub[krep];
			    maxdfs = xprune[krep];
	if (PROFILE) {
			    Gstat.procstat[pnum].pruned++;
	}
			} else {
			    fsupc = SUPER_FSUPC( xsup_end, supno[krep] );
			    xdfs = xlsub[fsupc] + krep-fsupc+1;
			    maxdfs = xlsub_end[fsupc];
	if (PROFILE) {
			    Gstat.procstat[pnum].unpruned++;
	}
			}

			do {
			    /*
			     * For each unmarked kchild of krep ...
			     */
			    while ( xdfs < maxdfs ) {

				kchild = lsub[xdfs];
				xdfs++;
				chmark = marker2[kchild];

				if ( chmark != jcol ) { /* Not reached yet */
				    marker2[kchild] = jcol;
				    chperm = perm_r[kchild];

				    if ( chperm == EMPTY ) {
					/* kchild is in L: place it in L[*,k]. */
					++no_lsub;
					col_lsub[nextl++] = kchild;
					if (chmark != jcolm1) samesuper = NO.ordinal();
				    } else {
					/* kchild is in U: chrep = its supernode
					 * representative. If its rep has
					 * been explored, update its repfnz[*].
					 */
					chrep = SUPER_REP( xsup_end, supno[chperm] );
					myfnz = repfnz[chrep];
					if ( myfnz != EMPTY ) { /* Visited before */
					    if ( myfnz > chperm )
						repfnz[chrep] = chperm;
					} else {
					    /* Continue dfs at super-rep of kchild */
					    xplore[krep] = xdfs;
					    xplore[m + krep] = maxdfs;
					    parent[chrep] = krep;
					    krep = chrep; /* Go deeper down G(L^t) */
					    repfnz[krep] = chperm;
					    if ( ispruned[krep] != 0 ) {
						if ( SINGLETON( xsup_end, xsup_end, supno[krep] ) )
						    xdfs = xlsub_end[krep];
						else xdfs = xlsub[krep];
						maxdfs = xprune[krep];
	if (PROFILE) {
						Gstat.procstat[pnum].pruned++;
	}
					    } else {
						fsupc = SUPER_FSUPC( xsup_end, supno[krep] );
						xdfs = xlsub[fsupc] + krep-fsupc+1;
						maxdfs = xlsub_end[fsupc];
	if (PROFILE) {
						Gstat.procstat[pnum].unpruned++;
	}
					    }
					}
				    } /* else */
				} /* if */
			    } /* while */

			    /* krow has no more unexplored nbrs:
			     *    place supernode-rep krep in postorder DFS,
			     *    backtrack dfs to its parent.
			     */
			    segrep[nseg[0]] = krep;
			    ++(nseg[0]);
	if ( DEBUGlevel>=3 ) {
	  if (jcol == BADCOL)
	    printf("(%d) pdgstrf_column_dfs[inner-dfs] new nseg %d, repfnz[krep=%d] %d\n",
		   pnum, nseg[0], krep, repfnz[krep]);
	}
			    kpar = parent[krep]; /* Pop from stack, mimic recursion */
			    if ( kpar == EMPTY ) break; /* dfs done */
			    krep = kpar;
			    xdfs = xplore[krep];
			    maxdfs = xplore[m + krep];
			} while ( kpar != EMPTY ); /* Do ... until empty stack */

		    } /* else myfnz ... */
		} /* if kperm >= fstcol ... */
	    } /* for each nonzero ... */

	if ( DEBUGlevel>=3 ) {
	  if (jcol == BADCOL)
	    printf("(%d) pdgstrf_column_dfs[2]: nextl %d, samesuper? %d\n",
		   pnum, nextl, samesuper);
	}

	    /* assert(no_lsub == nextl - no_usub);*/

	    /* ---------------------------------------------------------
	       Check to see if j belongs in the same supernode as j-1.
	       --------------------------------------------------------- */

	    /* Does it matter if jcol == 0? - xiaoye */
	    if ( samesuper == YES.ordinal() ) {
		nsuper = supno[jcolm1];
		jcolm1size = xlsub_end[jcolm1] - xlsub[jcolm1];
	if ( DEBUGlevel>=3 ) {
	  if (jcol == BADCOL)
	    printf("(%d) pdgstrf_column_dfs[YES] jcol-1 %d, jcolm1size %d, supno[%d] %d\n",
		   pnum, jcolm1, jcolm1size, jcolm1, nsuper);
	}
		if ( no_lsub != jcolm1size-1 )
		    samesuper = NO.ordinal();        /* Enforce T2 supernode */
		else {
		    /* Make sure the number of columns in a supernode does not
		       exceed threshold. */
		    fsupc = xsup[nsuper];
		    if ( jcol - fsupc >= maxsuper )
			samesuper = NO.ordinal();
		    else {
			/* start of a supernode in H (coarser partition) */
			if ( super_bnd[jcol] != 0 ) samesuper = NO.ordinal();
		    }
		}
	    }

	    /* If jcol starts a new supernode, allocate storage for
	     * the subscript set of both first and last column of
	     * a previous supernode. (first for num values, last for pruning)
	     */
	    if ( samesuper == NO.ordinal() ) { /* starts a new supernode */
	    int[] nsuper_ = new int[1];
		nsuper = NewNsuper(pnum, pxgstrf_shared, nsuper_);
		Glu.nsuper = nsuper_[0];
		xsup[nsuper] = jcol;

		/* Copy column jcol; also reserve space to store pruned graph */
		if ((mem_error = Glu_alloc(pnum, jcol, 2*no_lsub, LSUB, ito,
					  pxgstrf_shared)) != 0)
		    return mem_error;
		xlsub[jcol] = ito[0];
		lsub = Glu.lsub;
		for (ifrom = 0; ifrom < nextl; ++ifrom) {
		    krow = col_lsub[ifrom];
		    if ( perm_r[krow] == EMPTY ) /* Filter U-subscript */
			lsub[ito[0]++] = krow;
		}
		k = ito[0];
		xlsub_end[jcol] = k;

		/* Make a copy in case it is a singleton supernode */
		for (ifrom = xlsub[jcol]; ifrom < ito[0]; ++ifrom)
		    lsub[k++] = lsub[ifrom];

	    } else { /* Supernode of size > 1: overwrite column jcol-1 */
		k = xlsub_end[fsupc];
		xlsub[jcol] = k;
		xprune[fsupc] = k;
		for (ifrom = 0; ifrom < nextl; ++ifrom) {
		    krow = col_lsub[ifrom];
		    if ( perm_r[krow] == EMPTY ) /* Filter U-subscript */
			lsub[k++] = krow;
		}
		xlsub_end[jcol] = k;
	    }

	if ( DEBUGlevel>=3 ) {
	  if (jcol == BADCOL) {
	    printf("(%d) pdgstrf_column_dfs[3]: %d in prev s-node %d? %d\n",
		   pnum, jcol, fsupc, samesuper);
	    PrintInt10("lsub", xlsub_end[jcol]-xlsub[jcol], lsub, xlsub[jcol]);
	  }
	}

	    /* Tidy up the pointers before exit */
	    xprune[jcol] = k;     /* upper bound for pruning */
	    supno[jcol] = nsuper;
	    xsup_end[nsuper] = jcol + 1;

	    return 0;
	}

}

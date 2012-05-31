package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;

import static gov.lbl.superlu.Dlu.CHK_PRUNE;
import static gov.lbl.superlu.Dlu.printf;
import static gov.lbl.superlu.Dlu_slu_mt_util.EMPTY;
import static gov.lbl.superlu.Dlu_slu_mt_util.FALSE;
import static gov.lbl.superlu.Dlu_slu_mt_util.HICOL;
import static gov.lbl.superlu.Dlu_slu_mt_util.LOCOL;
import static gov.lbl.superlu.Dlu_slu_mt_util.SINGLETON;
import static gov.lbl.superlu.Dlu_slu_mt_util.TRUE;



public class Dlu_pxgstrf_pruneL {

	static
	void
	pxgstrf_pruneL(
		       final int  jcol,      /* current column */
		       final int  perm_r[],  /* row pivotings */
		       final int  pivrow,    /* pivot row of column jcol */
		       final int  nseg,      /* number of U-segments */
		       final int  segrep[],   /* in */
		       final int  repfnz[],   /* in */
		       int        xprune[],   /* modified */
		       int        ispruned[], /* modified */
		       GlobalLU_t Glu /* modified - global LU data structures */
		       )
	{
	/*
	 * -- SuperLU MT routine (version 1.0) --
	 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
	 * and Lawrence Berkeley National Lab.
	 * August 15, 1997
	 *
	 * Purpose
	 * =======
	 *   Reduces the L-structure of those supernodes whose L-structure
	 *   contains the current pivot row "pivrow".
	 *
	 */
	    int jsupno, irep, isupno, irep1, kmin = 0, kmax = 0, krow;
	    int i, ktemp;
	    int do_prune; /* logical variable */
	    int[]        xsup, xsup_end, supno;
	    int[]        lsub, xlsub, xlsub_end;

	    xsup       = Glu.xsup;
	    xsup_end   = Glu.xsup_end;
	    supno      = Glu.supno;
	    lsub       = Glu.lsub;
	    xlsub      = Glu.xlsub;
	    xlsub_end  = Glu.xlsub_end;

	    /*
	     * For each supernode-rep irep in U[*,j]
	     */
	    jsupno = supno[jcol];
	    for (i = 0; i < nseg; i++) {

		irep = segrep[i];
		irep1 = irep + 1;

		/* Don't prune with a zero U-segment */
	 	if ( repfnz[irep] == EMPTY ) continue;

	     	/* If a supernode overlaps with the next panel, then the U-segment
	   	 * is fragmented into two parts - irep and irep1. We should let
		 * pruning occur at the rep-column in irep1's supernode.
		 */
		isupno = supno[irep];
		if ( isupno == supno[irep1] ) continue;	/* Don't prune */

		/*
		 * If it is not pruned & it has a nonz in row L[pivrow,i]
		 */
		do_prune = FALSE;
		if ( isupno != jsupno ) {
		    if ( ispruned[irep] == 0 ) {
			kmin = SINGLETON( xsup_end, xsup_end, isupno ) ? xlsub_end[irep] : xlsub[irep];
			kmax = xprune[irep] - 1;
			for (krow = kmin; krow <= kmax; krow++)
			    if ( lsub[krow] == pivrow ) {
				do_prune = TRUE;
				break;
			    }
		    }

	    	    if ( do_prune != 0 ) {

		     	/* Do a quicksort-type partition */
		        while ( kmin <= kmax ) {
		    	    if ( perm_r[lsub[kmax]] == EMPTY )
				kmax--;
			    else if ( perm_r[lsub[kmin]] != EMPTY )
				kmin++;
			    else { /* kmin below pivrow, and kmax above pivrow:
			            * 	interchange the two subscripts
				    */
			        ktemp = lsub[kmin];
			        lsub[kmin] = lsub[kmax];
			        lsub[kmax] = ktemp;
			        kmin++;
			        kmax--;
			    }
		        } /* while */

		        xprune[irep] = kmin;	/* Pruning */
			ispruned[irep] = 1;

	if (CHK_PRUNE) {
	if (irep >= LOCOL && irep >= HICOL && jcol >= LOCOL && jcol <= HICOL)
	    printf("pxgstrf_pruneL() for irep %d using col %d: xprune %d - %d\n",
		   irep, jcol, xlsub[irep], kmin);
	}
		    } /* if do_prune */

		} /* if */

	    } /* for each U-segment... */
	}

}

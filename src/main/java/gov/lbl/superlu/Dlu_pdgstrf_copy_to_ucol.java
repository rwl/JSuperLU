package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_pdsp_defs.pxgstrf_shared_t;

import static gov.lbl.superlu.Dlu.DEBUG;
import static gov.lbl.superlu.Dlu_pmemory.Glu_alloc;
import static gov.lbl.superlu.Dlu_slu_mt_util.EMPTY;
import static gov.lbl.superlu.Dlu_slu_mt_util.MemType.UCOL;




public class Dlu_pdgstrf_copy_to_ucol {

	static
	int
	pdgstrf_copy_to_ucol(
			     final int  pnum,    /* process number */
			     final int  jcol,	 /* current column */
			     final int  nseg,	 /* number of U-segments */
			     final int  segrep[], /* in */
			     final int  repfnz[], /* in */
			     final int  perm_r[], /* in */
			     double	 dense[],  /* modified - reset to zero on exit */
			     pxgstrf_shared_t pxgstrf_shared /* modified */
			     )
	{
	/*
	 * -- SuperLU MT routine (version 2.0) --
	 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
	 * and Xerox Palo Alto Research Center.
	 * September 10, 2007
	 *
	 * Gather the nonzeros from SPA dense[*,jcol] into global ucol[*].
	 */
	    int ksub, krep, ksupno, i, k, kfnz, segsze;
	    int fsupc, isub, irow, jsupno, colsize;
	    int      nextu[], mem_error;
	    nextu = new int[1];
	    int[]    xsup, supno, lsub, xlsub, usub;
	    double   ucol[];
	    GlobalLU_t Glu = pxgstrf_shared.Glu; /* modified */

	    double zero = 0.0;

	    xsup    = Glu.xsup;
	    supno   = Glu.supno;
	    lsub    = Glu.lsub;
	    xlsub   = Glu.xlsub;
	    jsupno  = supno[jcol];

	    /* find the size of column jcol */
	    colsize = 0;
	    k = nseg - 1;
	    for (ksub = 0; ksub < nseg; ++ksub) {
		krep = segrep[k--];
		ksupno = supno[krep];

		if ( ksupno != jsupno ) { /* should go into ucol[] */
		    kfnz = repfnz[krep];
		    if ( kfnz != EMPTY )  /* nonzero U-segment */
			colsize += krep - kfnz + 1;;
		}
	    } /* for each segment... */

	    if ( (mem_error = Glu_alloc(pnum, jcol, colsize, UCOL, nextu,
					pxgstrf_shared)) != 0 )
		return mem_error;
	    Glu.xusub[jcol] = nextu[0];
	    ucol = Glu.ucol;
	    usub = Glu.usub;

	    /* Now, it does not have to be in topological order! */
	    k = nseg - 1;
	    for (ksub = 0; ksub < nseg; ++ksub) {

		krep = segrep[k--];
		ksupno = supno[krep];

		if ( ksupno != jsupno ) { /* should go into ucol[] */
		    kfnz = repfnz[krep];
		    if ( kfnz != EMPTY ) { /* nonzero U-segment */
		    	fsupc = xsup[ksupno];
		        isub = xlsub[fsupc] + kfnz - fsupc;
		        segsze = krep - kfnz + 1;
//	#pragma ivdep
			for (i = 0; i < segsze; i++) {
			    irow = lsub[isub];
			    usub[nextu[0]] = perm_r[irow];
			    ucol[nextu[0]] = dense[irow];
			    dense[irow] = zero;
	if (DEBUG) {
	if (jcol == EMPTY)
	    printf("(%d) pcopy_to_ucol[]: jcol %d, krep %d, irow %d, ucol %.10e\n",
		   ME, jcol, krep, irow, ucol[nextu[0]]);
	}
			    nextu[0]++;
			    isub++;
			}
		    }
		}

	    } /* for each segment... */

	    Glu.xusub_end[jcol] = nextu[0]; /* close U[*,jcol] */
	    return 0;
	}

}

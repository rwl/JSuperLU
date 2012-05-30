package lu.jsuper;

import org.netlib.blas.BLAS;

import lu.jsuper.Dlu_pdsp_defs.GlobalLU_t;
import lu.jsuper.Dlu_slu_mt_util.Gstat_t;

import static lu.jsuper.Dlu.USE_VENDOR_BLAS;

import static lu.jsuper.Dlu_dmyblas2.dlsolve;
import static lu.jsuper.Dlu_dmyblas2.dmatvec;


public class Dlu_pdgstrf_snode_bmod {

	static
	int
	pdgstrf_snode_bmod(
			   final int  pnum,   /* process number */
			   final int  jcol,   /* in - current column in the s-node */
			   final int  jsupno, /* in */
			   final int  fsupc,  /* in - first column in the s-node */
			   double     dense[], /* in */
			   double     tempv[], /* working array */
			   GlobalLU_t Glu,   /* modified */
			   Gstat_t Gstat     /* modified */
			   )
	{
	/*
	 * -- SuperLU MT routine (version 2.0) --
	 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
	 * and Xerox Palo Alto Research Center.
	 * September 10, 2007
	 *
	 * Performs numeric block updates within the relaxed supernode.
	 */

	    double      zero = 0.0;
	    double      one = 1.0;
	    double      none = -1.0;
	    int            incx = 1, incy = 1;
	    double         alpha = none, beta = one;

	    int            luptr, nsupc, nsupr, nrow;
	    int            isub, irow, i, iptr;
	    int            ufirst, nextlu;
	    double         lusup[];
	    int[]          lsub, xlsub, xlsub_end, xlusup, xlusup_end;
	    float flopcnt;

	    lsub       = Glu.lsub;
	    xlsub      = Glu.xlsub;
	    xlsub_end  = Glu.xlsub_end;
	    lusup      = Glu.lusup;
	    xlusup     = Glu.xlusup;
	    xlusup_end = Glu.xlusup_end;

	    nextlu = xlusup[jcol];

	    /*
	     *	Process the supernodal portion of L\U[*,j]
	     */
	    for (isub = xlsub[fsupc]; isub < xlsub_end[fsupc]; isub++) {
	  	irow = lsub[isub];
		lusup[nextlu] = dense[irow];
		dense[irow] = zero;
		++nextlu;
	    }

	    xlusup_end[jcol] = nextlu;

	    if ( fsupc < jcol ) {

		luptr = xlusup[fsupc];
		nsupr = xlsub_end[fsupc] - xlsub[fsupc];
		nsupc = jcol - fsupc;	/* Excluding jcol */
		ufirst = xlusup[jcol];	/* Points to the beginning of column
					   jcol in supernode L\U(jsupno). */
		nrow = nsupr - nsupc;

	        flopcnt = nsupc * (nsupc - 1) + 2 * nrow * nsupc;
		Gstat.procstat[pnum].fcops += flopcnt;

	/*	ops[TRSV] += nsupc * (nsupc - 1);
		ops[GEMV] += 2 * nrow * nsupc;    */

	if (USE_VENDOR_BLAS) {
		BLAS blas = BLAS.getInstance();

		blas.dtrsv( "L", "N", "U", &nsupc, &lusup[luptr], &nsupr,
		      &lusup[ufirst], &incx );
		blas.dgemv( "N", &nrow, &nsupc, &alpha, &lusup[luptr+nsupc], &nsupr,
			&lusup[ufirst], &incx, &beta, &lusup[ufirst+nsupc], &incy );
	} else {
		dlsolve ( nsupr, nsupc, &lusup[luptr], &lusup[ufirst] );
		dmatvec ( nsupr, nrow, nsupc, &lusup[luptr+nsupc],
			 &lusup[ufirst], &tempv[0] );

	        /* Scatter tempv[*] into lusup[*] */
		iptr = ufirst + nsupc;
		for (i = 0; i < nrow; i++) {
	            lusup[iptr++] -= tempv[i];
	            tempv[i] = zero;
		}
	}

	    }

	    return 0;
	}

}

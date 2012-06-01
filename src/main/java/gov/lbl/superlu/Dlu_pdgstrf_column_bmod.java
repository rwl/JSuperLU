package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_pdsp_defs.pxgstrf_shared_t;
import gov.lbl.superlu.Dlu_slu_mt_util.Gstat_t;

import static gov.lbl.superlu.Dlu_slu_mt_util.BADCOL;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_MAX;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPER_FSUPC;
import static gov.lbl.superlu.Dlu_slu_mt_util.MemType.LUSUP;

import static gov.lbl.superlu.Dlu_pmemory.Glu_alloc;

import static gov.lbl.superlu.Dlu.DEBUGlevel;
import static gov.lbl.superlu.Dlu.printf;
import static gov.lbl.superlu.Dlu.USE_VENDOR_BLAS;
import static gov.lbl.superlu.Dlu.DEBUG;
import static gov.lbl.superlu.Dlu.dtrsv;
import static gov.lbl.superlu.Dlu.dgemv;

import static gov.lbl.superlu.Dlu_dmyblas2.dlsolve;
import static gov.lbl.superlu.Dlu_dmyblas2.dmatvec;

import static gov.lbl.superlu.Dlu_pdutil.print_double_vec;


public class Dlu_pdgstrf_column_bmod {

	static
	int
	pdgstrf_column_bmod(
			    final int  pnum,   /* process number */
			    final int  jcol,   /* current column in the panel */
			    final int  fpanelc,/* first column in the panel */
			    final int  nseg,   /* number of s-nodes to update jcol */
			    int        segrep[],/* in */
			    int        repfnz[],/* in */
			    double     dense[], /* modified */
			    double     tempv[], /* working array */
			    pxgstrf_shared_t pxgstrf_shared, /* modified */
			    Gstat_t Gstat     /* modified */
			    )
	{
	/*
	 * -- SuperLU MT routine (version 2.0) --
	 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
	 * and Xerox Palo Alto Research Center.
	 * September 10, 2007
	 *
	 * Purpose:
	 * ========
	 *    Performs numeric block updates (sup-col) in topological order.
	 *    It features: col-col, 2cols-col, 3cols-col, and sup-col updates.
	 *    Special processing on the supernodal portion of L\U[*,j].
	 *
	 * Return value:
	 * =============
	 *      0 - successful return
	 *    > 0 - number of bytes allocated when run out of space
	 *
	 */
	    int         incx = 1, incy = 1;
	    double      alpha, beta;
	    GlobalLU_t Glu = pxgstrf_shared.Glu;   /* modified */

	    /* krep = representative of current k-th supernode
	     * fsupc = first supernodal column
	     * nsupc = no of columns in supernode
	     * nsupr = no of rows in supernode (used as leading dimension)
	     * luptr = location of supernodal LU-block in storage
	     * kfnz = first nonz in the k-th supernodal segment
	     * no_zeros = no of leading zeros in a supernodal U-segment
	     */
	    double	  ukj, ukj1, ukj2;
	    int lptr, kfnz, isub, irow, i, no_zeros;
	    int luptr, luptr1, luptr2;
	    int          fsupc, nsupc, nsupr, segsze;
	    int          nrow;	  /* No of rows in the matrix of matrix-vector */
	    int          jsupno, k, ksub, krep, krep_ind, ksupno;
	    int          ufirst, nextlu[];
	    nextlu = new int[1];
	    int          fst_col; /* First column within small LU update */
	    int          d_fsupc; /* Distance between the first column of the current
				     panel and the first column of the current snode.*/
	    int          xsup[], supno[];
	    int          lsub[], xlsub[], xlsub_end[];
	    double       lusup[];
	    int          xlusup[], xlusup_end[];
	    double       tempv1[];
	    int          tempv1_offset;
	    int          mem_error;
	    float flopcnt;

	    double      zero = 0.0;
	    double      one = 1.0;
	    double      none = -1.0;

	    xsup       = Glu.xsup;
	    supno      = Glu.supno;
	    lsub       = Glu.lsub;
	    xlsub      = Glu.xlsub;
	    xlsub_end  = Glu.xlsub_end;
	    lusup      = Glu.lusup;
	    xlusup     = Glu.xlusup;
	    xlusup_end = Glu.xlusup_end;
	    jsupno     = supno[jcol];

	    /*
	     * For each nonz supernode segment of U[*,j] in topological order
	     */
	    k = nseg - 1;
	    for (ksub = 0; ksub < nseg; ksub++) {

		krep = segrep[k];
		k--;
		ksupno = supno[krep];
	if ( DEBUGlevel>=2 ) {
	if (jcol==BADCOL)
	printf("(%d) pdgstrf_column_bmod[1]: %d, nseg %d, krep %d, jsupno %d, ksupno %d\n",
	       pnum, jcol, nseg, krep, jsupno, ksupno);
	}
		if ( jsupno != ksupno ) { /* Outside the rectangular supernode */

		    fsupc = xsup[ksupno];
		    fst_col = SUPERLU_MAX ( fsupc, fpanelc );

	  	    /* Distance from the current supernode to the current panel;
		       d_fsupc=0 if fsupc >= fpanelc. */
	  	    d_fsupc = fst_col - fsupc;

		    luptr = xlusup[fst_col] + d_fsupc;
		    lptr = xlsub[fsupc] + d_fsupc;
		    kfnz = repfnz[krep];
		    kfnz = SUPERLU_MAX ( kfnz, fpanelc );
		    segsze = krep - kfnz + 1;
		    nsupc = krep - fst_col + 1;
		    nsupr = xlsub_end[fsupc] - xlsub[fsupc]; /* Leading dimension */
		    nrow = nsupr - d_fsupc - nsupc;
		    krep_ind = lptr + nsupc - 1;

	            flopcnt = segsze * (segsze - 1) + 2 * nrow * segsze;
		    Gstat.procstat[pnum].fcops += flopcnt;

	if ( DEBUGlevel>=2 ) {
	if (jcol==BADCOL)
	printf("(%d) pdgstrf_column_bmod[2]: %d, krep %d, kfnz %d, segsze %d, d_fsupc %d," +
	"fsupc %d, nsupr %d, nsupc %d\n",
	       pnum, jcol, krep, kfnz, segsze, d_fsupc, fsupc, nsupr, nsupc);
	}


		    /*
		     * Case 1: Update U-segment of size 1 -- col-col update
		     */
		    if ( segsze == 1 ) {
		  	ukj = dense[lsub[krep_ind]];
			luptr += nsupr*(nsupc-1) + nsupc;

			for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
			    irow = lsub[i];
			    dense[irow] -=  ukj*lusup[luptr];
			    luptr++;
			}
		    } else if ( segsze <= 3 ) {
			ukj = dense[lsub[krep_ind]];
			luptr += nsupr*(nsupc-1) + nsupc-1;
			ukj1 = dense[lsub[krep_ind - 1]];
			luptr1 = luptr - nsupr;
			if ( segsze == 2 ) { /* Case 2: 2cols-col update */
			    ukj -= ukj1 * lusup[luptr1];
			    dense[lsub[krep_ind]] = ukj;
			    for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
			    	irow = lsub[i];
			    	luptr++;
			    	luptr1++;
			    	dense[irow] -= ( ukj*lusup[luptr]
						+ ukj1*lusup[luptr1] );
			    }
			} else { /* Case 3: 3cols-col update */
			    ukj2 = dense[lsub[krep_ind - 2]];
			    luptr2 = luptr1 - nsupr;
			    ukj1 -= ukj2 * lusup[luptr2-1];
			    ukj = ukj - ukj1*lusup[luptr1] - ukj2*lusup[luptr2];
			    dense[lsub[krep_ind]] = ukj;
			    dense[lsub[krep_ind-1]] = ukj1;
			    for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
			    	irow = lsub[i];
			    	luptr++;
			    	luptr1++;
				luptr2++;
			    	dense[irow] -= ( ukj*lusup[luptr]
				     + ukj1*lusup[luptr1] + ukj2*lusup[luptr2] );
			    }
			}


		    } else {
		  	/*
			 * Case: sup-col update
			 * Perform a triangular solve and block update,
			 * then scatter the result of sup-col update to dense
			 */
			no_zeros = kfnz - fst_col;

		        /* Copy U[*,j] segment from dense[*] to tempv[*] */
		        isub = lptr + no_zeros;
		        for (i = 0; i < segsze; i++) {
		  	    irow = lsub[isub];
			    tempv[i] = dense[irow];
			    ++isub;
		        }

		        /* Dense triangular solve -- start effective triangle */
			luptr += nsupr * no_zeros + no_zeros;
	if (USE_VENDOR_BLAS) {
			dtrsv( "L", "N", "U", segsze, lusup, luptr,
			       nsupr, tempv, 0, incx );

	 		luptr += segsze;  /* Dense matrix-vector */
			tempv1 = tempv;
			tempv1_offset = segsze;
			alpha = one;
			beta = zero;

			dgemv( "N", nrow, segsze, alpha, lusup, luptr,
			       nsupr, tempv, 0, incx, beta, tempv1, tempv1_offset, incy );
	} else {
			dlsolve ( nsupr, segsze, lusup, luptr, tempv, 0 );

	 		luptr += segsze;  /* Dense matrix-vector */
			tempv1 = tempv;
			tempv1_offset = segsze;
			dmatvec (nsupr, nrow , segsze, lusup, luptr, tempv, tempv1_offset, tempv1, 0);
	}
	                /* Scatter tempv[] into SPA dense[*] */
	                isub = lptr + no_zeros;
	                for (i = 0; i < segsze; i++) {
	                    irow = lsub[isub];
	                    dense[irow] = tempv[i]; /* Scatter */
	                    tempv[i] = zero;
	                    isub++;
	                }

			/* Scatter tempv1[] into SPA dense[*] */
			for (i = 0; i < nrow; i++) {
			    irow = lsub[isub];
	                    dense[irow] -= tempv1[i];
			    tempv1[i] = zero;
			    ++isub;
			}
		    } /* else segsze >= 4 */

		} /* if jsupno ... */

	    } /* for each segment... */


	    /* ------------------------------------------
	       Process the supernodal portion of L\U[*,j]
	       ------------------------------------------ */

	    fsupc = SUPER_FSUPC(xsup, jsupno);
	    nsupr = xlsub_end[fsupc] - xlsub[fsupc];
	    if ( (mem_error = Glu_alloc(pnum, jcol, nsupr, LUSUP, nextlu,
				       pxgstrf_shared)) != 0 )
		return mem_error;
	    xlusup[jcol] = nextlu[0];
	    lusup = Glu.lusup;

	    /* Gather the nonzeros from SPA dense[*,j] into L\U[*,j] */
	    for (isub = xlsub[fsupc]; isub < xlsub_end[fsupc]; ++isub) {
	  	irow = lsub[isub];
		lusup[nextlu[0]] = dense[irow];
		dense[irow] = zero;
	if (DEBUG) {
	if (jcol == -1)
	    printf("(%d) pdgstrf_column_bmod[lusup] jcol %d, irow %d, lusup %.10e\n",
		   pnum, jcol, irow, lusup[nextlu[0]]);
	}
		++nextlu[0];
	    }
	    xlusup_end[jcol] = nextlu[0]; /* close L\U[*,jcol] */

	if ( DEBUGlevel>=2 ) {
	if (jcol == -1) {
	    nrow = xlusup_end[jcol] - xlusup[jcol];
	    print_double_vec("before sup-col update", nrow, lsub, xlsub[fsupc],
			     lusup, xlusup[jcol]);
	}
	}

	    /*
	     * For more updates within the panel (also within the current supernode),
	     * should start from the first column of the panel, or the first column
	     * of the supernode, whichever is bigger. There are 2 cases:
	     *    (1) fsupc < fpanelc,  then fst_col := fpanelc
	     *    (2) fsupc >= fpanelc, then fst_col := fsupc
	     */
	    fst_col = SUPERLU_MAX ( fsupc, fpanelc );

	    if ( fst_col < jcol ) {

	  	/* distance between the current supernode and the current panel;
		   d_fsupc=0 if fsupc >= fpanelc. */
	  	d_fsupc = fst_col - fsupc;

		lptr = xlsub[fsupc] + d_fsupc;
		luptr = xlusup[fst_col] + d_fsupc;
		nsupr = xlsub_end[fsupc] - xlsub[fsupc]; /* Leading dimension */
		nsupc = jcol - fst_col;	/* Excluding jcol */
		nrow = nsupr - d_fsupc - nsupc;

		/* points to the beginning of jcol in supernode L\U[*,jsupno] */
		ufirst = xlusup[jcol] + d_fsupc;

	if ( DEBUGlevel>=2 ) {
	if (jcol==BADCOL)
	printf("(%d) pdgstrf_column_bmod[3] jcol %d, fsupc %d, nsupr %d, nsupc %d, nrow %d\n",
	       pnum, jcol, fsupc, nsupr, nsupc, nrow);
	}

	        flopcnt = nsupc * (nsupc - 1) + 2 * nrow * nsupc;
		Gstat.procstat[pnum].fcops += flopcnt;

	/*	ops[TRSV] += nsupc * (nsupc - 1);
		ops[GEMV] += 2 * nrow * nsupc;    */

	if (USE_VENDOR_BLAS) {
		alpha = none; beta = one; /* y := beta*y + alpha*A*x */
		dtrsv( "L", "N", "U", nsupc, lusup, luptr,
		       nsupr, lusup, ufirst, incx );
		dgemv( "N", nrow, nsupc, alpha, lusup, luptr+nsupc, nsupr,
		       lusup, ufirst, incx, beta, lusup, ufirst+nsupc, incy );
	} else {
		dlsolve ( nsupr, nsupc, lusup, luptr, lusup, ufirst );

		dmatvec ( nsupr, nrow, nsupc, lusup, luptr+nsupc,
			 lusup, ufirst, tempv, 0 );

	        /* Copy updates from tempv[*] into lusup[*] */
		isub = ufirst + nsupc;
		for (i = 0; i < nrow; i++) {
	            lusup[isub] -= tempv[i];
	            tempv[i] = 0.0;
		    ++isub;
		}
	}
	    } /* if fst_col < jcol ... */

	    return 0;
	}

}

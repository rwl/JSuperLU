/*! @file Dlu_dpanel_bmod.java
 * \brief Performs numeric block updates
 *
 * <pre>
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 * Copyright (c) 1994 by Xerox Corporation.  All rights reserved.
 *
 * THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
 * EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
 *
 * Permission is hereby granted to use or copy this program for any
 * purpose, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is
 * granted, provided the above notices are retained, and a notice that
 * the code was modified is included with the above copyright notice.
 * </pre>
 */
package lu.jsuper;

import org.netlib.blas.BLAS;

import lu.jsuper.Dlu_slu_ddefs.GlobalLU_t;
import lu.jsuper.Dlu_slu_util.SuperLUStat_t;
import lu.jsuper.Dlu_superlu_enum_consts.PhaseType;

import static lu.jsuper.Dlu_sp_ienv.sp_ienv;
import static lu.jsuper.Dlu_slu_util.SUPERLU_MAX;
import static lu.jsuper.Dlu_slu_util.SUPERLU_MIN;
import static lu.jsuper.Dlu_slu_util.EMPTY;
import static lu.jsuper.Dlu_util.USE_VENDOR_BLAS;

import static lu.jsuper.Dlu_dmyblas2.dlsolve;
import static lu.jsuper.Dlu_dmyblas2.dmatvec;


public class Dlu_dpanel_bmod {

    public static int first = 1, maxsuper, rowblk, colblk;

	/**! \brief
	 *
	 * <pre>
	 * Purpose
	 * =======
	 *
	 *    Performs numeric block updates (sup-panel) in topological order.
	 *    It features: col-col, 2cols-col, 3cols-col, and sup-col updates.
	 *    Special processing on the supernodal portion of L\U[*,j]
	 *
	 *    Before entering this routine, the original nonzeros in the panel
	 *    were already copied into the spa[m,w].
	 *
	 *    Updated/Output parameters-
	 *    dense[0:m-1,w]: L[*,j:j+w-1] and U[*,j:j+w-1] are returned
	 *    collectively in the m-by-w vector dense[*].
	 * </pre>
	 */
	public static void
	dpanel_bmod (
		    final int  m,          /* in - number of rows in the matrix */
		    final int  w,          /* in */
		    final int  jcol,       /* in */
		    final int  nseg,       /* in */
		    double     dense[],    /* out, of size n by w */
		    double     tempv[],    /* working array */
		    int        segrep[],   /* in */
		    int        repfnz[],   /* in, of size n by w */
		    GlobalLU_t Glu,        /* modified */
		    SuperLUStat_t stat     /* output */
		    )
	{
	    int          incx = 1, incy = 1;
	    double       alpha, beta;

	    int          k, ksub;
	    int          fsupc, nsupc, nsupr, nrow;
	    int          krep, krep_ind;
	    double       ukj, ukj1, ukj2;
	    int          luptr, luptr1, luptr2;
	    int          segsze;
	    int          block_nrow;  /* no of rows in a block row */
	    int          lptr;	      /* Points to the row subscripts of a supernode */
	    int          kfnz, irow, no_zeros;
	    int          isub, isub1, i;
	    int          jj;	      /* Index through each column in the panel */
	    int          xsup[], supno[];
	    int          lsub[], xlsub[];
	    double       lusup[];
	    int          xlusup[];
	    int          repfnz_col[]; /* repfnz[] for a column in the panel */
	    int          repfnz_col_offset;
	    double       dense_col[];  /* dense[] for a column in the panel */
	    int          dense_col_offset;
	    double       tempv1[];             /* Used in 1-D update */
	    int          tempv1_offset;
	    double       TriTmp[], MatvecTmp[]; /* used in 2-D update */
	    int          TriTmp_offset;
	    int          MatvecTmp_offset;
	    double       zero = 0.0;
	    double       one = 1.0;
	    int          ldaTmp;
	    int          r_ind, r_hi;
	    float        ops[] = stat.ops;

	    xsup    = Glu.xsup;
	    supno   = Glu.supno;
	    lsub    = Glu.lsub;
	    xlsub   = Glu.xlsub;
	    lusup   = Glu.lusup;
	    xlusup  = Glu.xlusup;

	    if ( first != 0 ) {
		maxsuper = SUPERLU_MAX( sp_ienv(3), sp_ienv(7) );
		rowblk   = sp_ienv(4);
		colblk   = sp_ienv(5);
		first = 0;
	    }
	    ldaTmp = maxsuper + rowblk;

	    /*
	     * For each nonz supernode segment of U[*,j] in topological order
	     */
	    k = nseg - 1;
	    for (ksub = 0; ksub < nseg; ksub++) { /* for each updating supernode */

		/* krep = representative of current k-th supernode
		 * fsupc = first supernodal column
		 * nsupc = no of columns in a supernode
		 * nsupr = no of rows in a supernode
		 */
	        krep = segrep[k--];
		fsupc = xsup[supno[krep]];
		nsupc = krep - fsupc + 1;
		nsupr = xlsub[fsupc+1] - xlsub[fsupc];
		nrow = nsupr - nsupc;
		lptr = xlsub[fsupc];
		krep_ind = lptr + nsupc - 1;

		repfnz_col = repfnz;
		repfnz_col_offset = 0;
		dense_col = dense;
		dense_col_offset = 0;

		if ( nsupc >= colblk && nrow > rowblk ) { /* 2-D block update */

		    TriTmp = tempv;
		    TriTmp_offset = 0;

		    /* Sequence through each column in panel -- triangular solves */
		    for (jj = jcol; jj < jcol + w; jj++,
			 repfnz_col_offset += m, dense_col_offset += m, TriTmp_offset += ldaTmp ) {

			kfnz = repfnz_col[repfnz_col_offset+krep];
			if ( kfnz == EMPTY ) continue;	/* Skip any zero segment */

			segsze = krep - kfnz + 1;
			luptr = xlusup[fsupc];

			ops[PhaseType.TRSV.ordinal()] += segsze * (segsze - 1);
			ops[PhaseType.GEMV.ordinal()] += 2 * nrow * segsze;

			/* Case 1: Update U-segment of size 1 -- col-col update */
			if ( segsze == 1 ) {
			    ukj = dense_col[dense_col_offset+lsub[krep_ind]];
			    luptr += nsupr*(nsupc-1) + nsupc;

			    for (i = lptr + nsupc; i < xlsub[fsupc+1]; i++) {
				irow = lsub[i];
				dense_col[dense_col_offset+irow] -= ukj * lusup[luptr];
				++luptr;
			    }

			} else if ( segsze <= 3 ) {
			    ukj = dense_col[dense_col_offset+lsub[krep_ind]];
			    ukj1 = dense_col[dense_col_offset+lsub[krep_ind - 1]];
			    luptr += nsupr*(nsupc-1) + nsupc-1;
			    luptr1 = luptr - nsupr;

			    if ( segsze == 2 ) {
				ukj -= ukj1 * lusup[luptr1];
				dense_col[dense_col_offset+lsub[krep_ind]] = ukj;
				for (i = lptr + nsupc; i < xlsub[fsupc+1]; ++i) {
				    irow = lsub[i];
				    luptr++; luptr1++;
				    dense_col[dense_col_offset+irow] -= (ukj*lusup[luptr]
							+ ukj1*lusup[luptr1]);
				}
			    } else {
				ukj2 = dense_col[dense_col_offset+lsub[krep_ind - 2]];
				luptr2 = luptr1 - nsupr;
				ukj1 -= ukj2 * lusup[luptr2-1];
				ukj = ukj - ukj1*lusup[luptr1] - ukj2*lusup[luptr2];
				dense_col[dense_col_offset+lsub[krep_ind]] = ukj;
				dense_col[dense_col_offset+lsub[krep_ind-1]] = ukj1;
				for (i = lptr + nsupc; i < xlsub[fsupc+1]; ++i) {
				    irow = lsub[i];
				    luptr++; luptr1++; luptr2++;
				    dense_col[dense_col_offset+irow] -= ( ukj*lusup[luptr]
	                             + ukj1*lusup[luptr1] + ukj2*lusup[luptr2] );
				}
			    }

			} else  {	/* segsze >= 4 */

			    /* Copy U[*,j] segment from dense[*] to TriTmp[*], which
			       holds the result of triangular solves.    */
			    no_zeros = kfnz - fsupc;
			    isub = lptr + no_zeros;
			    for (i = 0; i < segsze; ++i) {
				irow = lsub[isub];
				TriTmp[TriTmp_offset+i] = dense_col[dense_col_offset+irow]; /* Gather */
				++isub;
			    }

			    /* start effective triangle */
			    luptr += nsupr * no_zeros + no_zeros;

			    if (USE_VENDOR_BLAS) {
			    BLAS blas = BLAS.getInstance();
			    blas.dtrsv( "L", "N", "U", segsze, lusup[luptr],
				   nsupr, TriTmp, TriTmp_offset, incx );
			    } else {
			    dlsolve ( nsupr, segsze, lusup[luptr], TriTmp, TriTmp_offset );
			    }


			} /* else ... */

		    }  /* for jj ... end tri-solves */

		    /* Block row updates; push all the way into dense[*] block */
		    for ( r_ind = 0; r_ind < nrow; r_ind += rowblk ) {

			r_hi = SUPERLU_MIN(nrow, r_ind + rowblk);
			block_nrow = SUPERLU_MIN(rowblk, r_hi - r_ind);
			luptr = xlusup[fsupc] + nsupc + r_ind;
			isub1 = lptr + nsupc + r_ind;

			repfnz_col = repfnz;
			repfnz_col_offset = 0;
			TriTmp = tempv;
			TriTmp_offset = 0;
			dense_col = dense;
			dense_col_offset = 0;

			/* Sequence through each column in panel -- matrix-vector */
			for (jj = jcol; jj < jcol + w; jj++,
			     repfnz_col_offset += m, dense_col_offset += m, TriTmp_offset += ldaTmp) {

			    kfnz = repfnz_col[repfnz_col_offset+krep];
			    if ( kfnz == EMPTY ) continue; /* Skip any zero segment */

			    segsze = krep - kfnz + 1;
			    if ( segsze <= 3 ) continue;   /* skip unrolled cases */

			    /* Perform a block update, and scatter the result of
			       matrix-vector to dense[].		 */
			    no_zeros = kfnz - fsupc;
			    luptr1 = luptr + nsupr * no_zeros;
			    MatvecTmp = TriTmp;
			    MatvecTmp_offset = TriTmp_offset+maxsuper;

			    if (USE_VENDOR_BLAS) {
			    alpha = one;
	            beta = zero;
	            BLAS blas = BLAS.getInstance();
			    blas.dgemv("N", &block_nrow, segsze, alpha, lusup[luptr1],
				   nsupr, TriTmp, TriTmp_offset, incx, beta, MatvecTmp,
				   MatvecTmp_offset, incy);
			    } else {
			    dmatvec(nsupr, block_nrow, segsze, &lusup[luptr1],
				   TriTmp, TriTmp_offset, MatvecTmp, MatvecTmp_offset);
			    }

			    /* Scatter MatvecTmp[*] into SPA dense[*] temporarily
			     * such that MatvecTmp[*] can be re-used for the
			     * the next blok row update. dense[] will be copied into
			     * global store after the whole panel has been finished.
			     */
			    isub = isub1;
			    for (i = 0; i < block_nrow; i++) {
				irow = lsub[isub];
				dense_col[dense_col_offset+irow] -= MatvecTmp[MatvecTmp_offset+i];
				MatvecTmp[MatvecTmp_offset+i] = zero;
				++isub;
			    }

			} /* for jj ... */

		    } /* for each block row ... */

		    /* Scatter the triangular solves into SPA dense[*] */
		    repfnz_col = repfnz;
		    repfnz_col_offset = 0;
		    TriTmp = tempv;
		    TriTmp_offset = 0;
		    dense_col = dense;
		    dense_col_offset = 0;

		    for (jj = jcol; jj < jcol + w; jj++,
			 repfnz_col_offset += m, dense_col_offset += m, TriTmp_offset += ldaTmp) {
			kfnz = repfnz_col[repfnz_col_offset+krep];
			if ( kfnz == EMPTY ) continue; /* Skip any zero segment */

			segsze = krep - kfnz + 1;
			if ( segsze <= 3 ) continue; /* skip unrolled cases */

			no_zeros = kfnz - fsupc;
			isub = lptr + no_zeros;
			for (i = 0; i < segsze; i++) {
			    irow = lsub[isub];
			    dense_col[dense_col_offset+irow] = TriTmp[TriTmp_offset+i];
			    TriTmp[TriTmp_offset+i] = zero;
			    ++isub;
			}

		    } /* for jj ... */

		} else { /* 1-D block modification */


		    /* Sequence through each column in the panel */
		    for (jj = jcol; jj < jcol + w; jj++,
			 repfnz_col_offset += m, dense_col_offset += m) {

			kfnz = repfnz_col[repfnz_col_offset+krep];
			if ( kfnz == EMPTY ) continue;	/* Skip any zero segment */

			segsze = krep - kfnz + 1;
			luptr = xlusup[fsupc];

			ops[PhaseType.TRSV.ordinal()] += segsze * (segsze - 1);
			ops[PhaseType.GEMV.ordinal()] += 2 * nrow * segsze;

			/* Case 1: Update U-segment of size 1 -- col-col update */
			if ( segsze == 1 ) {
			    ukj = dense_col[dense_col_offset+lsub[krep_ind]];
			    luptr += nsupr*(nsupc-1) + nsupc;

			    for (i = lptr + nsupc; i < xlsub[fsupc+1]; i++) {
				irow = lsub[i];
				dense_col[dense_col_offset+irow] -= ukj * lusup[luptr];
				++luptr;
			    }

			} else if ( segsze <= 3 ) {
			    ukj = dense_col[dense_col_offset+lsub[krep_ind]];
			    luptr += nsupr*(nsupc-1) + nsupc-1;
			    ukj1 = dense_col[dense_col_offset+lsub[krep_ind - 1]];
			    luptr1 = luptr - nsupr;

			    if ( segsze == 2 ) {
				ukj -= ukj1 * lusup[luptr1];
				dense_col[dense_col_offset+lsub[krep_ind]] = ukj;
				for (i = lptr + nsupc; i < xlsub[fsupc+1]; ++i) {
				    irow = lsub[i];
				    ++luptr;  ++luptr1;
				    dense_col[dense_col_offset+irow] -= (ukj*lusup[luptr]
							+ ukj1*lusup[luptr1]);
				}
			    } else {
				ukj2 = dense_col[dense_col_offset+lsub[krep_ind - 2]];
				luptr2 = luptr1 - nsupr;
				ukj1 -= ukj2 * lusup[luptr2-1];
				ukj = ukj - ukj1*lusup[luptr1] - ukj2*lusup[luptr2];
				dense_col[dense_col_offset+lsub[krep_ind]] = ukj;
				dense_col[dense_col_offset+lsub[krep_ind-1]] = ukj1;
				for (i = lptr + nsupc; i < xlsub[fsupc+1]; ++i) {
				    irow = lsub[i];
				    ++luptr; ++luptr1; ++luptr2;
				    dense_col[dense_col_offset+irow] -= ( ukj*lusup[luptr]
	                             + ukj1*lusup[luptr1] + ukj2*lusup[luptr2] );
				}
			    }

			} else  { /* segsze >= 4 */
			    /*
			     * Perform a triangular solve and block update,
			     * then scatter the result of sup-col update to dense[].
			     */
			    no_zeros = kfnz - fsupc;

			    /* Copy U[*,j] segment from dense[*] to tempv[*]:
			     *    The result of triangular solve is in tempv[*];
			     *    The result of matrix vector update is in dense_col[*]
			     */
			    isub = lptr + no_zeros;
			    for (i = 0; i < segsze; ++i) {
				irow = lsub[isub];
				tempv[i] = dense_col[dense_col_offset+irow]; /* Gather */
				++isub;
			    }

			    /* start effective triangle */
			    luptr += nsupr * no_zeros + no_zeros;

			    if (USE_VENDOR_BLAS) {
			    BLAS blas = BLAS.getInstance();
			    blas.dtrsv( "L", "N", "U", segsze, lusup[luptr],
				   nsupr, tempv, incx );


			    luptr += segsze;	/* Dense matrix-vector */
			    tempv1 = tempv;
			    tempv1_offset = segsze;
                alpha = one;
                beta = zero;

			    blas.dgemv( "N", nrow, segsze, alpha, lusup[luptr],
				   nsupr, tempv, incx, beta, tempv1, tempv1_offset, incy );
			    }
			    else
			    {
			    dlsolve ( nsupr, segsze, lusup[luptr], tempv );

			    luptr += segsze;        /* Dense matrix-vector */
			    tempv1 = tempv;
			    tempv1_offset = segsze;
			    dmatvec (nsupr, nrow, segsze, &lusup[luptr], tempv, tempv1, tempv1_offset);
		        }

			    /* Scatter tempv[*] into SPA dense[*] temporarily, such
			     * that tempv[*] can be used for the triangular solve of
			     * the next column of the panel. They will be copied into
			     * ucol[*] after the whole panel has been finished.
			     */
			    isub = lptr + no_zeros;
			    for (i = 0; i < segsze; i++) {
				irow = lsub[isub];
				dense_col[dense_col_offset+irow] = tempv[i];
				tempv[i] = zero;
				isub++;
			    }

			    /* Scatter the update from tempv1[*] into SPA dense[*] */
			    /* Start dense rectangular L */
			    for (i = 0; i < nrow; i++) {
				irow = lsub[isub];
				dense_col[dense_col_offset+irow] -= tempv1[tempv1_offset+i];
				tempv1[tempv1_offset+i] = zero;
				++isub;
			    }

			} /* else segsze>=4 ... */

		    } /* for each column in the panel... */

		} /* else 1-D update ... */

	    } /* for each updating supernode ... */

	}

}

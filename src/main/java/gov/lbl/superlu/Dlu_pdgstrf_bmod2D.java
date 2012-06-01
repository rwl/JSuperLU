package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_slu_mt_util.Gstat_t;

import static gov.lbl.superlu.Dlu_sp_ienv.sp_ienv;

import static gov.lbl.superlu.Dlu_slu_mt_util.EMPTY;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.FLOAT;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_MIN;

import static gov.lbl.superlu.Dlu_superlu_timer.SuperLU_timer_;

import static gov.lbl.superlu.Dlu.TIMING;
import static gov.lbl.superlu.Dlu.SCATTER_FOUND;
import static gov.lbl.superlu.Dlu.USE_VENDOR_BLAS;
import static gov.lbl.superlu.Dlu.dtrsv;
import static gov.lbl.superlu.Dlu.dgemv;

import static gov.lbl.superlu.Dlu_dmyblas2.dlsolve;
import static gov.lbl.superlu.Dlu_dmyblas2.dmatvec;


public class Dlu_pdgstrf_bmod2D {

    static int first = 1, maxsuper, rowblk;

	static
	void
	pdgstrf_bmod2D(
		       final int pnum,   /* process number */
		       final int m,      /* number of columns in the matrix */
		       final int w,      /* current panel width */
		       final int jcol,   /* leading column of the current panel */
		       final int fsupc,  /* leading column of the updating supernode */
		       final int krep,   /* last column of the updating supernode */
		       final int nsupc,  /* number of columns in the updating s-node */
		       int nsupr,        /* number of rows in the updating s-node */
		       int nrow,         /* number of rows below the diagonal block of
					    the updating supernode */
		       int repfnz[],      /* in */
		       int panel_lsub[],  /* modified */
		       int w_lsub_end[],  /* modified */
		       int spa_marker[],  /* modified; size n-by-w */
		       double dense[],    /* modified */
		       double tempv[],    /* working array - zeros on entry/exit */
		       GlobalLU_t Glu,  /* modified */
		       Gstat_t Gstat    /* modified */
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
	 *
	 *    Performs numeric 2-D block updates (sup-panel) in topological order.
	 *    Results are returned in SPA dense[*,w].
	 *
	 */
	    int          incx = 1, incy = 1;
	    double      alpha, beta;

	    double      zero = 0.0;
	    double      one = 1.0;

	    double       ukj, ukj1, ukj2;
	    int          luptr, luptr1, luptr2;
	    int          segsze;
	    int          block_nrow;  /* no of rows in a block row */
	    int lptr;	      /* points to the row subscripts of a supernode */
	    int          kfnz, irow, no_zeros;
	    int isub, isub1, i;
	    int jj;	      /* index through each column in the panel */
	    int          krep_ind;
	    int          repfnz_col[]; /* repfnz[] for a column in the panel */
	    int          col_marker[]; /* each column of the spa_marker[*,w] */
	    int          col_lsub[];   /* each column of the panel_lsub[*,w] */
	    double       dense_col[];  /* dense[] for a column in the panel */
	    double       TriTmp[], MatvecTmp[];
	    int ldaTmp;
	    int r_ind, r_hi;
	    int          lsub[], xlsub_end[];
	    double       lusup[];
	    int          xlusup[];
	    float flopcnt;

	    double utime[] = Gstat.utime;
	    double f_time;

	    if ( first != 0 ) {
		maxsuper = sp_ienv(3);
		rowblk   = sp_ienv(4);
		first = 0;
	    }
	    ldaTmp = maxsuper + rowblk;

	    lsub      = Glu.lsub;
	    xlsub_end = Glu.xlsub_end;
	    lusup     = Glu.lusup;
	    xlusup    = Glu.xlusup;
	    lptr      = Glu.xlsub[fsupc];
	    krep_ind  = lptr + nsupc - 1;
	    repfnz_col= repfnz;
	    dense_col = dense;
	    TriTmp    = tempv;
	    col_marker= spa_marker;
	    col_lsub  = panel_lsub;


	    /* ---------------------------------------------------------------
	     * Sequence through each column in the panel -- triangular solves.
	     * The results of the triangular solves of all columns in the
	     * panel are temporaroly stored in TriTemp array.
	     * For the unrolled small supernodes of size <= 3, we also perform
	     * matrix-vector updates from below the diagonal block.
	     * ---------------------------------------------------------------
	     */
	    for (jj = jcol; jj < jcol + w; ++jj, col_marker += m, col_lsub += m,
		 repfnz_col += m, dense_col += m, TriTmp += ldaTmp ) {

		kfnz = repfnz_col[krep];
		if ( kfnz == EMPTY ) continue;	/* Skip any zero segment */

		segsze = krep - kfnz + 1;
		luptr = xlusup[fsupc];

	        flopcnt = segsze * (segsze - 1) + 2 * nrow * segsze;
		Gstat.procstat[pnum].fcops += flopcnt;

	/*	ops[TRSV] += segsze * (segsze - 1);
		ops[GEMV] += 2 * nrow * segsze;        */

	if (TIMING) {
		f_time = SuperLU_timer_();
	}

		/* Case 1: Update U-segment of size 1 -- col-col update */
		if ( segsze == 1 ) {
		    ukj = dense_col[lsub[krep_ind]];
		    luptr += nsupr*(nsupc-1) + nsupc;
		    for (i = lptr + nsupc; i < xlsub_end[fsupc]; i++) {
			irow = lsub[i];
	                dense_col[irow] -= ukj * lusup[luptr];
			++luptr;
	if (SCATTER_FOUND) {
			if ( col_marker[irow] != jj ) {
			    col_marker[irow] = jj;
			    col_lsub[w_lsub_end[jj-jcol]++] = irow;
			}
	}
		    }
	if (TIMING) {
		    utime[FLOAT.ordinal()] += SuperLU_timer_() - f_time;
	}
		} else if ( segsze <= 3 ) {
		    ukj = dense_col[lsub[krep_ind]];
		    ukj1 = dense_col[lsub[krep_ind - 1]];
		    luptr += nsupr*(nsupc-1) + nsupc-1;
		    luptr1 = luptr - nsupr;
		    if ( segsze == 2 ) {
	                ukj -= ukj1 * lusup[luptr1];
			dense_col[lsub[krep_ind]] = ukj;
			for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
			    irow = lsub[i];
			    luptr++; luptr1++;
	                    dense_col[irow] -= (ukj * lusup[luptr]
	                                                + ukj1 * lusup[luptr1]);
	if (SCATTER_FOUND) {
			    if ( col_marker[irow] != jj ) {
				col_marker[irow] = jj;
				col_lsub[w_lsub_end[jj-jcol]++] = irow;
			    }
	}
			}
	if (TIMING) {
			utime[FLOAT.ordinal()] += SuperLU_timer_() - f_time;
	}
		    } else {
			ukj2 = dense_col[lsub[krep_ind - 2]];
			luptr2 = luptr1 - nsupr;
	                ukj1 -= ukj2 * lusup[luptr2-1];
	                ukj = ukj - ukj1*lusup[luptr1] - ukj2*lusup[luptr2];
			dense_col[lsub[krep_ind]] = ukj;
			dense_col[lsub[krep_ind-1]] = ukj1;
			for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
			    irow = lsub[i];
			    luptr++; luptr1++; luptr2++;
	                    dense_col[irow] -= (ukj * lusup[luptr]
	                             + ukj1*lusup[luptr1] + ukj2*lusup[luptr2]);
	if (SCATTER_FOUND) {
			    if ( col_marker[irow] != jj ) {
				col_marker[irow] = jj;
				col_lsub[w_lsub_end[jj-jcol]++] = irow;
			    }
	}
			}
		    }
	if (TIMING) {
		    utime[FLOAT.ordinal()] += SuperLU_timer_() - f_time;
	}
		} else  { /* segsze >= 4 */
		    /* Copy A[*,j] segment from dense[*] to TriTmp[*], which
		       holds the result of triangular solve.    */
		    no_zeros = kfnz - fsupc;
		    isub = lptr + no_zeros;
		    for (i = 0; i < segsze; ++i) {
			irow = lsub[isub];
			TriTmp[i] = dense_col[irow]; /* Gather */
			++isub;
		    }

		    /* start effective triangle */
		    luptr += nsupr * no_zeros + no_zeros;

	if (TIMING) {
		    f_time = SuperLU_timer_();
	}

	if (USE_VENDOR_BLAS) {
		    dtrsv( "L", "N", "U", segsze, lusup, luptr,
			   nsupr, TriTmp, 0, incx );
	} else {
		    dlsolve ( nsupr, segsze, lusup, luptr, TriTmp, 0 );
	}

	if (TIMING) {
		    utime[FLOAT.ordinal()] += SuperLU_timer_() - f_time;
	}
		} /* else ... */

	    }  /* for jj ... end tri-solves */

	    /* --------------------------------------------------------
	     * Perform block row updates from below the diagonal block.
	     * Push each block all the way into SPA dense[*].
	     * --------------------------------------------------------
	     */
	    for ( r_ind = 0; r_ind < nrow; r_ind += rowblk ) {

		r_hi = SUPERLU_MIN(nrow, r_ind + rowblk);
		block_nrow = SUPERLU_MIN(rowblk, r_hi - r_ind);
		luptr = xlusup[fsupc] + nsupc + r_ind;
		isub1 = lptr + nsupc + r_ind;

		repfnz_col = repfnz;
		TriTmp = tempv;
		dense_col = dense;
		col_marker= spa_marker;
		col_lsub  = panel_lsub;

		/* Sequence through each column in the panel -- matrix-vector */
		for (jj = jcol; jj < jcol + w; ++jj, col_marker += m, col_lsub += m,
		     repfnz_col += m, dense_col += m, TriTmp += ldaTmp) {

		    kfnz = repfnz_col[krep];
		    if ( kfnz == EMPTY ) continue; /* skip any zero segment */

		    segsze = krep - kfnz + 1;
		    if ( segsze <= 3 ) continue;   /* skip unrolled cases */

		    /* Perform a block update, and scatter the result of
		       matrix-vector into SPA dense[*].		 */
		    no_zeros = kfnz - fsupc;
		    luptr1 = luptr + nsupr * no_zeros;
		    MatvecTmp = &TriTmp[maxsuper];

	if (TIMING) {
		    f_time = SuperLU_timer_();
	}

	if (USE_VENDOR_BLAS) {
	            alpha = one;
	            beta = zero;
		    dgemv( "N", block_nrow, segsze, alpha, lusup, luptr1,
			   nsupr, TriTmp, 0, incx, beta, MatvecTmp, 0, incy );
	} else {
		    dmatvec(nsupr, block_nrow, segsze, lusup, luptr1,
			    TriTmp, 0, MatvecTmp);
	}

	if (TIMING) {
		    utime[FLOAT.ordinal()] += SuperLU_timer_() - f_time;
	}

		    /* Scatter MatvecTmp[*] into SPA dense[*] temporarily,
		     * such that MatvecTmp[*] can be re-used for the
		     * the next block row update. dense[] will be copied into
		     * global store after the whole panel has been finished.
		     */
		    isub = isub1;
		    for (i = 0; i < block_nrow; i++) {
			irow = lsub[isub];
	                dense_col[irow] -= MatvecTmp[i]; /* Scatter-add */
	if (SCATTER_FOUND) {
			if ( col_marker[irow] != jj ) {
			    col_marker[irow] = jj;
			    col_lsub[w_lsub_end[jj-jcol]++] = irow;
			}
	}
			MatvecTmp[i] = zero;
			++isub;
		    }

		} /* for jj ... */

	    } /* for each block row ... */


	    /* ------------------------------------------------
	       Scatter the triangular solves into SPA dense[*].
	       ------------------------------------------------ */
	    repfnz_col = repfnz;
	    TriTmp = tempv;
	    dense_col = dense;

	    for (jj = 0; jj < w; ++jj, repfnz_col += m, dense_col += m,
		 TriTmp += ldaTmp) {
		kfnz = repfnz_col[krep];
		if ( kfnz == EMPTY ) continue; /* skip any zero segment */

		segsze = krep - kfnz + 1;
		if ( segsze <= 3 ) continue; /* skip unrolled cases */

		no_zeros = kfnz - fsupc;
		isub = lptr + no_zeros;
		for (i = 0; i < segsze; i++) {
		    irow = lsub[isub];
		    dense_col[irow] = TriTmp[i]; /* Scatter */
		    TriTmp[i] = zero;
		    ++isub;
		}
	    } /* for jj ... */

	}
}

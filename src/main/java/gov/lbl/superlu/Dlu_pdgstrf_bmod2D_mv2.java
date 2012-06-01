package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_slu_mt_util.Gstat_t;

import static gov.lbl.superlu.Dlu_sp_ienv.sp_ienv;

import static gov.lbl.superlu.Dlu_slu_mt_util.EMPTY;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.FLOAT;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_MIN;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_MAX;

import static gov.lbl.superlu.Dlu.TIMING;
import static gov.lbl.superlu.Dlu.SCATTER_FOUND;
import static gov.lbl.superlu.Dlu.USE_VENDOR_BLAS;
import static gov.lbl.superlu.Dlu.dtrsv;
import static gov.lbl.superlu.Dlu.dgemv;

import static gov.lbl.superlu.Dlu_superlu_timer.SuperLU_timer_;

import static gov.lbl.superlu.Dlu_dmyblas2.dlsolve;


public class Dlu_pdgstrf_bmod2D_mv2 {

    static int first = 1, maxsuper, rowblk;

	static
	void
	pdgstrf_bmod2D_mv2(
			   final int pnum, /* process number */
			   final int n,    /* number of rows in the matrix */
			   final int w,    /* current panel width */
			   final int jcol, /* leading column of the current panel */
			   final int fsupc,/* leading column of the updating supernode */
			   final int krep, /* last column of the updating s-node */
			   final int nsupc,/* number of columns in the updating s-node */
			   int nsupr, /* number of rows in the updating s-node */
			   int nrow,  /* number of rows below the diagonal block of
					 the updating supernode */
			   int repfnz[],     /* in */
			   int panel_lsub[], /* modified */
			   int w_lsub_end[], /* modified */
			   int spa_marker[], /* modified; size n-by-w */
			   double dense[],   /* modified */
			   double tempv[],   /* working array - zeros on entry/exit */
			   GlobalLU_t Glu, /* modified */
			   Gstat_t Gstat   /* modified */
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

	    double      zero = 0.0;
	    double      one = 1.0;

	    int          incx = 1, incy = 1;
	    double       alpha = one, beta = zero;

	    double       ukj, ukj1, ukj2;
	    int          luptr, luptr1, luptr2;
	    int          segsze;
	    int          block_nrow;  /* no of rows in a block row */
	    int lptr; /* points to the row subscripts of a supernode */
	    int          kfnz, irow, no_zeros;
	    int isub, isub1, i, j;
	    int jj;	      /* index through each column in the panel */
	    int          krep_ind;
	    int          repfnz_col[]; /* repfnz[] for a column in the panel */
	    int          col_marker[]; /* each column of the spa_marker[*,w] */
	    int          col_lsub[];   /* each column of the panel_lsub[*,w] */
	    double       dense_col[];  /* dense[] for a column in the panel */
	    double       TriTmp[];
	    int ldaTmp;
	    int r_ind, r_hi;
	    int twocols;
	    int          kfnz2[], jj2[]; /* detect two identical columns */
	    kfnz2 = new int[2];
	    jj2 = new int[2];
	    double       tri[], matvec[];
	    tri = new double[2];
	    matvec = new double[2];
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


	    /* ---------------------------------------------------------------
	     * Sequence through each column in the panel -- triangular solves.
	     * The results of the triangular solves of all columns in the
	     * panel are temporaroly stored in TriTemp[*] array.
	     * For the unrolled small supernodes of size <= 3, we also perform
	     * matrix-vector updates from below the diagonal block.
	     * ---------------------------------------------------------------
	     */
	    repfnz_col= repfnz;
	    dense_col = dense;
	    TriTmp    = tempv;
	    col_marker= spa_marker;
	    col_lsub  = panel_lsub;

	    for (jj = jcol; jj < jcol + w; ++jj, col_marker += n, col_lsub += n,
		 repfnz_col += n, dense_col += n, TriTmp += ldaTmp ) {

		kfnz = repfnz_col[krep];
		if ( kfnz == EMPTY ) continue;	/* Skip any zero segment */

		segsze = krep - kfnz + 1;
		luptr = xlusup[fsupc];

	        flopcnt = segsze * (segsze - 1) + 2 * nrow * segsze;
		Gstat.procstat[pnum].fcops += flopcnt;

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
		luptr1 = xlusup[fsupc] + nsupc + r_ind;
		isub1 = lptr + nsupc + r_ind;
		repfnz_col = repfnz;
		twocols = 0;

		/* Sequence through each column in the panel -- matrix-vector */
		for (jj = jcol; jj < jcol + w; ++jj, repfnz_col += n) {

		    kfnz = repfnz_col[krep];
		    if ( kfnz == EMPTY ) continue; /* skip zero segment */
		    segsze = krep - kfnz + 1;
		    if ( segsze <= 3 ) continue;   /* skip unrolled cases */

		    /* Now segsze >= 4 ... */

		    if ( twocols == 1 ) { /* got two columns */
			jj2[1] = jj;
			twocols = 0;
			for (j = 0; j < 2; ++j) {
			    i = n * (jj2[j] - jcol);
			    kfnz2[j] = repfnz[i + krep];
			    tri[j] = tempv + ldaTmp * (jj2[j] - jcol);
			    matvec[j] = tri[j] + maxsuper;
			}

			if ( kfnz2[0] < kfnz2[1] ) { /* First column is bigger */
			    no_zeros = kfnz2[0] - fsupc;
			    segsze = kfnz2[1] - kfnz2[0];
			    luptr = luptr1 + nsupr * no_zeros;
	if (USE_VENDOR_BLAS) {
			    dgemv( "N", block_nrow, segsze, alpha, lusup, luptr,
				   nsupr, tri, 0, incx, beta, matvec, 0, incy );
	} else {
			    dmatvec (nsupr, block_nrow, segsze, lusup, luptr,
				     tri, 0, matvec, 0);
	}
			} else if ( kfnz2[0] > kfnz2[1] ) {
			    no_zeros = kfnz2[1] - fsupc;
			    segsze = kfnz2[0] - kfnz2[1];
			    luptr = luptr1 + nsupr * no_zeros;
	if (USE_VENDOR_BLAS) {
			    dgemv( "N", block_nrow, segsze, alpha, lusup, luptr,
				   nsupr, tri, 1, incx, beta, matvec, 1, incy );
	} else {
			    dmatvec (nsupr, block_nrow, segsze, lusup, luptr,
				     tri, 1, matvec, 1);
	}
			}

			/* Do matrix-vector multiply with two destinations */
			kfnz = SUPERLU_MAX( kfnz2[0], kfnz2[1] );
			no_zeros = kfnz - fsupc;
			segsze = krep - kfnz + 1;
			luptr = luptr1 + nsupr * no_zeros;

			dmatvec2 (nsupr, block_nrow, segsze, &lusup[luptr],
				  &tri[0][kfnz-kfnz2[0]], &tri[1][kfnz-kfnz2[1]],
				  matvec[0], matvec[1]);

	if (TIMING) {
			utime[FLOAT.ordinal()] += SuperLU_timer_() - f_time;
	}
			/* end for two destination update */
		    } else { /* wait for a second column */
			jj2[0] = jj;
			twocols = 1;
		    }

		} /* for jj ... */

		if ( twocols == 1 ) { /* one more column left */
		    i = jj2[0] - jcol;
		    tri[0] = tempv + ldaTmp * i;
		    matvec[0] = tri[0] + maxsuper;
		    kfnz = repfnz[i*n + krep];
		    no_zeros = kfnz - fsupc;
		    segsze = krep - kfnz + 1;
		    luptr = luptr1 + nsupr * no_zeros;

	if (USE_VENDOR_BLAS) {
		    dgemv( "N", block_nrow, segsze, alpha, lusup, luptr,
			   nsupr, tri, 0, incx, beta, matvec, 0, incy );
	} else {
		    dmatvec(nsupr, block_nrow, segsze, lusup, luptr,
			    tri, 0, matvec, 0);
	}
		} /* if twocols == 1 */

		/* Scatter matvec[*] into SPA dense[*]. */
		repfnz_col = repfnz;
		dense_col = dense;
		col_marker = spa_marker;
		col_lsub = panel_lsub;
		matvec[0] = tempv + maxsuper;
		for (jj = jcol; jj < jcol + w; ++jj, repfnz_col += n, dense_col += n,
		     col_marker += n, col_lsub += n, matvec[0] += ldaTmp) {
		    kfnz = repfnz_col[krep];
		    if ( kfnz == EMPTY ) continue; /* skip zero segment */
		    segsze = krep - kfnz + 1;
		    if ( segsze <= 3 ) continue;   /* skip unrolled cases */

		    isub = isub1;
		    for (i = 0; i < block_nrow; ++i) {
			irow = lsub[isub];
	                dense_col[irow] -= matvec[0][i]; /* Scatter-add */
	if (SCATTER_FOUND) {
			if ( col_marker[irow] != jj ) {
			    col_marker[irow] = jj;
			    col_lsub[w_lsub_end[jj-jcol]++] = irow;
			}
	}
			matvec[0][i] = zero;
			++isub;
		    }
		} /* for jj ... */

	    } /* for each block row ... */


	    /* ------------------------------------------------
	       Scatter the triangular solves into SPA dense[*].
	       ------------------------------------------------ */
	    repfnz_col = repfnz;
	    dense_col = dense;
	    TriTmp = tempv;
	    for (jj = 0; jj < w; ++jj, repfnz_col += n,
		 dense_col += n, TriTmp += ldaTmp) {
		kfnz = repfnz_col[krep];
		if ( kfnz == EMPTY ) continue; /* skip any zero segment */
		segsze = krep - kfnz + 1;
		if ( segsze <= 3 ) continue;   /* skip unrolled cases */

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

package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_slu_mt_util.Gstat_t;

import static gov.lbl.superlu.Dlu_sp_ienv.sp_ienv;

import static gov.lbl.superlu.Dlu_slu_mt_util.EMPTY;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.FLOAT;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_MIN;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_MAX;
import static gov.lbl.superlu.Dlu_slu_mt_util.BADPAN;
import static gov.lbl.superlu.Dlu_slu_mt_util.BADREP;
import static gov.lbl.superlu.Dlu_slu_mt_util.BADCOL;

import static gov.lbl.superlu.Dlu_util.PrintInt10;

import static gov.lbl.superlu.Dlu.TIMING;
import static gov.lbl.superlu.Dlu.SCATTER_FOUND;
import static gov.lbl.superlu.Dlu.USE_VENDOR_BLAS;
import static gov.lbl.superlu.Dlu.dtrsv;
import static gov.lbl.superlu.Dlu.dgemv;
import static gov.lbl.superlu.Dlu.DEBUG;
import static gov.lbl.superlu.Dlu.printf;

import static gov.lbl.superlu.Dlu_superlu_timer.SuperLU_timer_;

import static gov.lbl.superlu.Dlu_dmyblas2.dlsolve;
import static gov.lbl.superlu.Dlu_dmyblas2.dmatvec;
import static gov.lbl.superlu.Dlu_dmyblas2.dmatvec2;


public class Dlu_pdgstrf_bmod1D_mv2 {

	static
	void
	pdgstrf_bmod1D_mv2(
			   final int pnum, /* process number */
			   final int n,    /* number of rows in the matrix */
			   final int w,    /* current panel width */
			   final int jcol, /* leading column of the current panel */
			   final int fsupc,/* leading column of the updating s-node */
			   final int krep, /* last column of the updating s-node */
			   final int nsupc,/* number of columns in the updating s-node */
			   int nsupr, /* number of rows in the updating supernode */
			   int nrow,  /* number of rows below the diagonal block of
					 the updating supernode */
			   int repfnz[],    /* in */
			   int panel_lsub[],/* modified */
			   int w_lsub_end[],/* modified */
			   int spa_marker[],/* modified; size n-by-w */
			   double dense[],  /* modified */
			   double tempv[],  /* working array - zeros on entry/exit */
			   GlobalLU_t Glu,/* modified */
			   Gstat_t Gstat  /* modified */
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
	 *    Performs numeric block updates (sup-panel) in topological order.
	 *    It features: col-col, 2cols-col, 3cols-col, and sup-col updates.
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
	    int lptr; /* start of row subscripts of the updating supernode */
	    int i, j, kfnz, krep_ind, isub, irow, no_zeros, twocols;
	    int jj;	       /* index through each column in the panel */
	    int      kfnz2[], jj2[]; /* detect two identical columns */
	    kfnz2 = new int[2];
	    jj2 = new int[2];
	    int  repfnz_col[], repfnz_col1[]; /* repfnz[] for a column in the panel */
	    double dense_col[], dense_col1[];  /* dense[] for a column in the panel */
	    double[] tri[], matvec[];
	    tri = new double[2][];
	    matvec = new double[2][];
	    int[] matvec_offset = {0, 0};
	    int  col_marker[], col_marker1[]; /* each column of the spa_marker[*,w] */
	    int  col_lsub[], col_lsub1[];   /* each column of the panel_lsub[*,w] */
	    int          lsub[], xlsub_end[];
	    double	lusup[];
	    int          xlusup[];
	    float flopcnt;

	    double utime[] = Gstat.utime;
	    double f_time = 0;

	    lsub      = Glu.lsub;
	    xlsub_end = Glu.xlsub_end;
	    lusup     = Glu.lusup;
	    xlusup    = Glu.xlusup;
	    lptr      = Glu.xlsub[fsupc];
	    krep_ind  = lptr + nsupc - 1;
	    twocols = 0;
	    tri[0] = tempv;
	    tri[1] = tempv;
	    int[] tri_offset = {0, n};

	if (DEBUG) {
	if (jcol == BADPAN && krep == BADREP) {
	    printf("(%d) dbmod1D[1] jcol %d, fsupc %d, krep %d, nsupc %d, nsupr %d, nrow %d\n",
		   pnum, jcol, fsupc, krep, nsupc, nsupr, nrow);
	    PrintInt10("lsub[xlsub[2774]", nsupr, lsub, lptr);
	}
	}

	    /* -----------------------------------------------
	     * Sequence through each column in the panel ...
	     * ----------------------------------------------- */
	    repfnz_col= repfnz;
	    dense_col = dense;
	    col_marker= spa_marker;
	    col_lsub  = panel_lsub;
	    int repfnz_col_offset = 0, dense_col_offset = 0;
	    int col_marker_offset = 0, col_lsub_offset = 0;

	    for (jj = jcol; jj < jcol + w; ++jj, col_marker_offset += n, col_lsub_offset += n,
		 repfnz_col_offset += n, dense_col_offset += n) {

		kfnz = repfnz_col[repfnz_col_offset+krep];
		if ( kfnz == EMPTY ) continue;	/* skip any zero segment */

		segsze = krep - kfnz + 1;
		luptr = xlusup[fsupc];

	        flopcnt = segsze * (segsze - 1) + 2 * nrow * segsze;
		Gstat.procstat[pnum].fcops += flopcnt;

		/* Case 1: Update U-segment of size 1 -- col-col update */
		if ( segsze == 1 ) {
	if (TIMING) {
		    f_time = SuperLU_timer_();
	}
		    ukj = dense_col[dense_col_offset + lsub[krep_ind]];
		    luptr += nsupr*(nsupc-1) + nsupc;
	if (DEBUG) {
	if (krep == BADCOL && jj == -1) {
	    printf("(%d) dbmod1D[segsze=1]: k %d, j %d, ukj %.10e\n",
		   pnum, lsub[krep_ind], jj, ukj);
	    PrintInt10("segsze=1", nsupr, lsub, lptr);
	}
	}
		    for (i = lptr + nsupc; i < xlsub_end[fsupc]; i++) {
			irow = lsub[i];
	                dense_col[dense_col_offset + irow] -= ukj * lusup[luptr];
			++luptr;
	if (SCATTER_FOUND) {
			if ( col_marker[col_marker_offset + irow] != jj ) {
			    col_marker[col_marker_offset + irow] = jj;
			    col_lsub[col_lsub_offset + w_lsub_end[jj-jcol]++] = irow;
			}
	}
		    }
	if (TIMING) {
		    utime[FLOAT.ordinal()] += SuperLU_timer_() - f_time;
	}

		} else if ( segsze <= 3 ) {

	if (TIMING) {
		    f_time = SuperLU_timer_();
	}
		    ukj = dense_col[dense_col_offset + lsub[krep_ind]];
		    luptr += nsupr*(nsupc-1) + nsupc-1;
		    ukj1 = dense_col[dense_col_offset + lsub[krep_ind - 1]];
		    luptr1 = luptr - nsupr;
		    if ( segsze == 2 ) {
	                ukj -= ukj1 * lusup[luptr1];
			dense_col[dense_col_offset + lsub[krep_ind]] = ukj;
	/*#pragma ivdep*/
			for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
			    irow = lsub[i];
			    ++luptr;  ++luptr1;
	                    dense_col[dense_col_offset + irow] -= (ukj * lusup[luptr]
	                                                + ukj1 * lusup[luptr1]);
	if (SCATTER_FOUND) {
			    if ( col_marker[col_marker_offset + irow] != jj ) {
				col_marker[col_marker_offset + irow] = jj;
				col_lsub[col_lsub_offset + w_lsub_end[jj-jcol]++] = irow;
			    }
	}
			}
		    } else {
			ukj2 = dense_col[dense_col_offset + lsub[krep_ind - 2]];
			luptr2 = luptr1 - nsupr;
	                ukj1 -= ukj2 * lusup[luptr2-1];
	                ukj = ukj - ukj1*lusup[luptr1] - ukj2*lusup[luptr2];
			dense_col[dense_col_offset + lsub[krep_ind]] = ukj;
			dense_col[dense_col_offset + lsub[krep_ind-1]] = ukj1;
			for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
			    irow = lsub[i];
			    ++luptr; ++luptr1; ++luptr2;
	                    dense_col[dense_col_offset + irow] -= (ukj * lusup[luptr]
	                             + ukj1*lusup[luptr1] + ukj2*lusup[luptr2]);
	if (SCATTER_FOUND) {
			    if ( col_marker[col_marker_offset + irow] != jj ) {
				col_marker[col_marker_offset + irow] = jj;
				col_lsub[col_lsub_offset + w_lsub_end[jj-jcol]++] = irow;
			    }
	}
			}
		    }
	if (TIMING) {
		    utime[FLOAT.ordinal()] += SuperLU_timer_() - f_time;
	}
		} else { /* segsze >= 4 */
		    if ( twocols == 1 ) {
			jj2[1] = jj; /* got two columns */
			twocols = 0;

			for (j = 0; j < 2; ++j) { /* Do two tri-solves */
			    i = n * (jj2[j] - jcol);
			    repfnz_col1 = repfnz;
			    int repfnz_col1_offset = i;
			    dense_col1  = dense;
			    int dense_col1_offset = i;
			    kfnz2[j] = repfnz_col1[repfnz_col1_offset+krep];
			    no_zeros = kfnz2[j] - fsupc;
			    segsze = krep - kfnz2[j] + 1;
			    matvec[j] = tri[j];
			    matvec_offset[j] = tri_offset[j] + segsze;

			    /* Gather U[*,j] segment from dense[*] to tri[*]. */
			    isub = lptr + no_zeros;
			    for (i = 0; i < segsze; ++i) {
				irow = lsub[isub];
				tri[j][tri_offset[j] + i] = dense_col1[dense_col1_offset+irow]; /* Gather */
				++isub;
			    }

	if (TIMING) {
			    f_time = SuperLU_timer_();
	}
			    /* start effective triangle */
			    luptr = xlusup[fsupc] + nsupr * no_zeros + no_zeros;

	if (USE_VENDOR_BLAS) {
			    dtrsv( "L", "N", "U", segsze, lusup, luptr,
				   nsupr, tri[j], tri_offset[j], incx );
	} else {
			    dlsolve ( nsupr, segsze, lusup, luptr, tri[j], tri_offset[j] );
	}

	if (TIMING) {
			    utime[FLOAT.ordinal()] += SuperLU_timer_() - f_time;
	}
			} /* end for j ... two tri-solves */

	if (TIMING) {
			f_time = SuperLU_timer_();
	}

			if ( kfnz2[0] < kfnz2[1] ) { /* First column is bigger */
			    no_zeros = kfnz2[0] - fsupc;
			    segsze = kfnz2[1] - kfnz2[0];
			    luptr = xlusup[fsupc] + nsupr * no_zeros + nsupc;
	if (USE_VENDOR_BLAS) {
			    dgemv( "N", nrow, segsze, alpha, lusup, luptr,
				   nsupr, tri[0], tri_offset[0], incx, beta, matvec[0], matvec_offset[0], incy );
	} else {
			    dmatvec (nsupr, nrow, segsze, lusup, luptr,
				     tri[0], tri_offset[0], matvec[0], matvec_offset[0]);
	}

			} else if ( kfnz2[0] > kfnz2[1] ) {
			    no_zeros = kfnz2[1] - fsupc;
			    segsze = kfnz2[0] - kfnz2[1];
			    luptr = xlusup[fsupc] + nsupr * no_zeros + nsupc;
	if (USE_VENDOR_BLAS) {
			    dgemv( "N", nrow, segsze, alpha, lusup, luptr,
				   nsupr, tri[1], tri_offset[1], incx, beta, matvec[1], matvec_offset[1], incy );
	} else {
			    dmatvec (nsupr, nrow, segsze, lusup, luptr,
				     tri[1], tri_offset[1], matvec[1], matvec_offset[1]);
	}
			}

			/* Do matrix-vector multiply with two destinations */
			kfnz = SUPERLU_MAX( kfnz2[0], kfnz2[1] );
			no_zeros = kfnz - fsupc;
			segsze = krep - kfnz + 1;
			luptr = xlusup[fsupc] + nsupr * no_zeros + nsupc;

			dmatvec2(nsupr, nrow, segsze, lusup, luptr,
				 tri[0], tri_offset[0]+kfnz-kfnz2[0],
				 tri[1], tri_offset[1]+kfnz-kfnz2[1],
				 matvec[0], matvec_offset[0], matvec[1], matvec_offset[1]);

	if (TIMING) {
			utime[FLOAT.ordinal()] += SuperLU_timer_() - f_time;
	}

			for (j = 0; j < 2; ++j) {
			    i = n * (jj2[j] - jcol);
			    dense_col1  = dense;
			    int dense_col1_offset = i;
			    col_marker1 = spa_marker;
			    int col_marker1_offset = i;
			    col_lsub1   = panel_lsub;
			    int col_lsub1_offset = i;
			    no_zeros = kfnz2[j] - fsupc;
			    segsze = krep - kfnz2[j] + 1;

			    /* Scatter tri[*] into SPA dense[*]. */
			    isub = lptr + no_zeros;
			    for (i = 0; i < segsze; i++) {
				irow = lsub[isub];
				dense_col1[dense_col1_offset+irow] = tri[j][tri_offset[j] + i]; /* Scatter */
				tri[j][tri_offset[j] + i] = zero;
				++isub;
	if (DEBUG) {
		if (jj == -1 && krep == 3423)
		    printf("(%d) dbmod1D[scatter] jj %d, dense_col[%d] %e\n",
			   pnum, jj, irow, dense_col[dense_col_offset + irow]);
	}
			    }

			    /* Scatter matvec[*] into SPA dense[*]. */
	/*#pragma ivdep*/
			    for (i = 0; i < nrow; i++) {
				irow = lsub[isub];
	            dense_col1[dense_col1_offset+irow] -= matvec[j][matvec_offset[j] + i]; /* Scatter-add */
	if (SCATTER_FOUND) {
				if ( col_marker1[col_marker1_offset+irow] != jj2[j] ) {
				    col_marker1[col_marker1_offset+irow] = jj2[j];
				    col_lsub1[col_lsub1_offset + w_lsub_end[jj2[j]-jcol]++] = irow;
				}
	}
				matvec[j][matvec_offset[j] + i] = zero;
				++isub;
			    }

			} /* end for two destination update */

		    } else { /* wait for a second column */
			jj2[0] = jj;
			twocols = 1;
		    }
		} /* else segsze >= 4 */

	    } /* for jj ... */


	    if ( twocols == 1 ) { /* one more column left */
		i = n * (jj2[0] - jcol);
		repfnz_col1 = repfnz;
		int repfnz_col1_offset = i;
		dense_col1  = dense;
		int dense_col1_offset = i;
		col_marker1 = spa_marker;
		int col_marker1_offset = i;
		col_lsub1   = panel_lsub;
		int col_lsub1_offset = i;
		kfnz = repfnz_col1[repfnz_col1_offset+krep];
		no_zeros = kfnz - fsupc;
		segsze = krep - kfnz + 1;

		/* Gather U[*,j] segment from dense[*] to tri[*]. */
		isub = lptr + no_zeros;
		for (i = 0; i < segsze; ++i) {
		    irow = lsub[isub];
		    tri[0][tri_offset[0] + i] = dense_col1[dense_col1_offset+irow]; /* Gather */
		    ++isub;
		}

	if (TIMING) {
		f_time = SuperLU_timer_();
	}
		/* start effective triangle */
		luptr = xlusup[fsupc] + nsupr * no_zeros + no_zeros;
	if (USE_VENDOR_BLAS) {
		dtrsv( "L", "N", "U", segsze, lusup, luptr,
		       nsupr, tri[0], tri_offset[0], incx );
	} else {
		dlsolve ( nsupr, segsze, lusup, luptr, tri[0], tri_offset[0] );
	}

		luptr += segsze;	/* Dense matrix-vector */
		matvec[0] = tri[0];
		matvec_offset[0] = tri_offset[0] + segsze;

	if (USE_VENDOR_BLAS) {
		dgemv( "N", nrow, segsze, alpha, lusup, luptr,
		       nsupr, tri[0], tri_offset[0], incx, beta, matvec[0], matvec_offset[0], incy );
	} else {
		dmatvec (nsupr, nrow, segsze, lusup, luptr, tri[0], tri_offset[0], matvec[0], matvec_offset[0]);
	}
	if (TIMING) {
		utime[FLOAT.ordinal()] += SuperLU_timer_() - f_time;
	}

		/* Scatter tri[*] into SPA dense[*]. */
		isub = lptr + no_zeros;
		for (i = 0; i < segsze; i++) {
		    irow = lsub[isub];
		    dense_col1[dense_col1_offset+irow] = tri[0][tri_offset[0] + i]; /* Scatter */
		    tri[0][tri_offset[0] + i] = zero;
		    ++isub;
	if (DEBUG) {
		if (jj == -1 && krep == 3423)
		    printf("(%d) dbmod1D[scatter] jj %d, dense_col[%d] %e\n",
			   pnum, jj, irow, dense_col[dense_col_offset + irow]);
	}
		}

		/* Scatter matvec[*] into SPA dense[*]. */
		for (i = 0; i < nrow; i++) {
		    irow = lsub[isub];
	        dense_col1[dense_col1_offset+irow] -= matvec[0][matvec_offset[0] + i]; /* Scatter-add */
	if (SCATTER_FOUND) {
		    if ( col_marker1[col_marker1_offset+irow] != jj2[0] ) {
			col_marker1[col_marker1_offset+irow] = jj2[0];
			col_lsub1[col_lsub1_offset + w_lsub_end[jj2[0]-jcol]++] = irow;
		    }
	}
		    matvec[0][matvec_offset[0] + i] = zero;
		    ++isub;
		}

	    } /* if twocols == 1 */

	}

}

package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_slu_mt_util.Gstat_t;

import static gov.lbl.superlu.Dlu_sp_ienv.sp_ienv;

import static gov.lbl.superlu.Dlu_slu_mt_util.EMPTY;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.FLOAT;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_MIN;
import static gov.lbl.superlu.Dlu_slu_mt_util.BADPAN;
import static gov.lbl.superlu.Dlu_slu_mt_util.BADREP;
import static gov.lbl.superlu.Dlu_slu_mt_util.BADCOL;

import static gov.lbl.superlu.Dlu_util.PrintInt10;

import static gov.lbl.superlu.Dlu_superlu_timer.SuperLU_timer_;

import static gov.lbl.superlu.Dlu.TIMING;
import static gov.lbl.superlu.Dlu.SCATTER_FOUND;
import static gov.lbl.superlu.Dlu.USE_VENDOR_BLAS;
import static gov.lbl.superlu.Dlu.dtrsv;
import static gov.lbl.superlu.Dlu.dgemv;
import static gov.lbl.superlu.Dlu.DEBUGlevel;
import static gov.lbl.superlu.Dlu.printf;

import static gov.lbl.superlu.Dlu_dmyblas2.dlsolve;
import static gov.lbl.superlu.Dlu_dmyblas2.dmatvec;


public class Dlu_pdgstrf_bmod1D {

	static
	void
	pdgstrf_bmod1D(
		       final int pnum,  /* process number */
		       final int m,     /* number of rows in the matrix */
		       final int w,     /* current panel width */
		       final int jcol,  /* leading column of the current panel */
		       final int fsupc, /* leading column of the updating supernode */
		       final int krep,  /* last column of the updating supernode */
		       final int nsupc, /* number of columns in the updating s-node */
		       int nsupr, /* number of rows in the updating supernode */
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
	 * Lawrence Berkeley National Lab,  Univ. of California Berkeley,
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
	    int          incx = 1, incy = 1;
	    double       alpha, beta;

	    double       ukj, ukj1, ukj2;
	    int          luptr, luptr1, luptr2;
	    int          segsze;
	    int lptr; /* start of row subscripts of the updating supernode */
	    int i, krep_ind, kfnz, isub, irow, no_zeros;
	    int jj;	      /* index through each column in the panel */
	    int          repfnz_col[]; /* repfnz[] for a column in the panel */
	    double       dense_col[];  /* dense[] for a column in the panel */
	    double       tempv1[];     /* used to store matrix-vector result */
	    int tempv1_offset;
	    int          col_marker[]; /* each column of the spa_marker[*,w] */
	    int          col_lsub[];   /* each column of the panel_lsub[*,w] */
	    int          lsub[], xlsub_end[];
	    double       lusup[];
	    int          xlusup[];
	    float flopcnt;

	    double      zero = 0.0;
	    double      one = 1.0;

	    double utime[] = Gstat.utime;
	    double f_time = 0;

	    lsub      = Glu.lsub;
	    xlsub_end = Glu.xlsub_end;
	    lusup     = Glu.lusup;
	    xlusup    = Glu.xlusup;
	    lptr      = Glu.xlsub[fsupc];
	    krep_ind  = lptr + nsupc - 1;

	    /* Pointers to each column of the w-wide arrays. */
	    repfnz_col= repfnz;
	    dense_col = dense;
	    col_marker= spa_marker;
	    col_lsub  = panel_lsub;
	    int repfnz_col_offset = 0, dense_col_offset = 0;
	    int col_marker_offset = 0, col_lsub_offset = 0;

	if ( DEBUGlevel>=2 ) {
	if (jcol == BADPAN && krep == BADREP) {
	    printf("(%d) pdgstrf_bmod1D[1] jcol %d, fsupc %d, krep %d, nsupc %d, nsupr %d, nrow %d\n",
		   pnum, jcol, fsupc, krep, nsupc, nsupr, nrow);
	    PrintInt10("lsub[xlsub[2774]]", nsupr, lsub, lptr);
	}
	}

	    /*
	     * Sequence through each column in the panel ...
	     */
	    for (jj = jcol; jj < jcol + w; ++jj, col_marker_offset += m, col_lsub_offset += m,
		 repfnz_col_offset += m, dense_col_offset += m) {

		kfnz = repfnz_col[repfnz_col_offset+krep];
		if ( kfnz == EMPTY ) continue;	/* Skip any zero segment */

		segsze = krep - kfnz + 1;
		luptr = xlusup[fsupc];

		/* Calculate flops: tri-solve + mat-vector */
	        flopcnt = segsze * (segsze - 1) + 2 * nrow * segsze;
		Gstat.procstat[pnum].fcops += flopcnt;

		/* Case 1: Update U-segment of size 1 -- col-col update */
		if ( segsze == 1 ) {
	if (TIMING) {
		    f_time = SuperLU_timer_();
	}
		    ukj = dense_col[dense_col_offset+lsub[krep_ind]];
		    luptr += nsupr*(nsupc-1) + nsupc;
	if ( DEBUGlevel>=2 ) {
	if (krep == BADCOL && jj == -1) {
	    printf("(%d) pdgstrf_bmod1D[segsze=1]: k %d, j %d, ukj %.10e\n",
		   pnum, lsub[krep_ind], jj, ukj);
	    PrintInt10("segsze=1", nsupr, lsub, lptr);
	}
	}
		    for (i = lptr + nsupc; i < xlsub_end[fsupc]; i++) {
			irow = lsub[i];
	                        dense_col[dense_col_offset+irow] -= ukj * lusup[luptr];
			++luptr;
	if (SCATTER_FOUND) {
			if ( col_marker[col_marker_offset+irow] != jj ) {
			    col_marker[col_marker_offset+irow] = jj;
			    col_lsub[col_lsub_offset+w_lsub_end[jj-jcol]++] = irow;
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
		    ukj = dense_col[dense_col_offset+lsub[krep_ind]];
		    luptr += nsupr*(nsupc-1) + nsupc-1;
		    ukj1 = dense_col[dense_col_offset+lsub[krep_ind - 1]];
		    luptr1 = luptr - nsupr;
		    if ( segsze == 2 ) {
	                ukj -= ukj1 * lusup[luptr1];
			dense_col[dense_col_offset+lsub[krep_ind]] = ukj;
			for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
			    irow = lsub[i];
			    ++luptr;  ++luptr1;
	                            dense_col[dense_col_offset+irow] -= (ukj * lusup[luptr]
	                                                + ukj1 * lusup[luptr1]);
	if (SCATTER_FOUND) {
			    if ( col_marker[col_marker_offset+irow] != jj ) {
				col_marker[col_marker_offset+irow] = jj;
				col_lsub[col_lsub_offset+w_lsub_end[jj-jcol]++] = irow;
			    }
	}
			}
		    } else {
			ukj2 = dense_col[dense_col_offset+lsub[krep_ind - 2]];
			luptr2 = luptr1 - nsupr;
	                ukj1 -= ukj2 * lusup[luptr2-1];
	                ukj = ukj - ukj1*lusup[luptr1] - ukj2*lusup[luptr2];
			dense_col[dense_col_offset+lsub[krep_ind]] = ukj;
			dense_col[dense_col_offset+lsub[krep_ind-1]] = ukj1;
			for (i = lptr + nsupc; i < xlsub_end[fsupc]; ++i) {
			    irow = lsub[i];
			    ++luptr; ++luptr1; ++luptr2;
	                    dense_col[dense_col_offset+irow] -= (ukj * lusup[luptr]
	                             + ukj1*lusup[luptr1] + ukj2*lusup[luptr2]);
	if (SCATTER_FOUND) {
			    if ( col_marker[col_marker_offset+irow] != jj ) {
				col_marker[col_marker_offset+irow] = jj;
				col_lsub[col_lsub_offset+w_lsub_end[jj-jcol]++] = irow;
			    }
	}
			}
		    }
	if (TIMING) {
		    utime[FLOAT.ordinal()] += SuperLU_timer_() - f_time;
	}
		} else { /* segsze >= 4 */
		    /*
		     * Perform a triangular solve and matrix-vector update,
		     * then scatter the result of sup-col update to dense[*].
		     */
		    no_zeros = kfnz - fsupc;

		    /* Gather U[*,j] segment from dense[*] to tempv[*]:
		     *   The result of triangular solve is in tempv[*];
		     *   The result of matrix vector update is in dense_col[*]
		     */
		    isub = lptr + no_zeros;
	/*#pragma ivdep*/
		    for (i = 0; i < segsze; ++i) {
			irow = lsub[isub];
			tempv[i] = dense_col[dense_col_offset+irow]; /* Gather */
			++isub;
		    }

		    /* start effective triangle */
		    luptr += nsupr * no_zeros + no_zeros;
	if (TIMING) {
		    f_time = SuperLU_timer_();
	}

	if (USE_VENDOR_BLAS) {
		    dtrsv( "L", "N", "U", segsze, lusup, luptr,
			   nsupr, tempv, 0, incx );

		    luptr += segsze;	/* Dense matrix-vector */
		    tempv1 = tempv;
		    tempv1_offset = segsze;

	        alpha = one;
	        beta = zero;

		    dgemv( "N", nrow, segsze, alpha, lusup, luptr,
			   nsupr, tempv, 0, incx, beta, tempv1, tempv1_offset, incy );
	} else {
		    dlsolve ( nsupr, segsze, lusup, luptr, tempv, 0 );

		    luptr += segsze;        /* Dense matrix-vector */
		    tempv1 = tempv;
		    tempv1_offset = segsze;
		    dmatvec (nsupr, nrow, segsze, lusup, luptr, tempv, 0, tempv1, tempv1_offset);
	}

	if (TIMING) {
		    utime[FLOAT.ordinal()] += SuperLU_timer_() - f_time;
	}

		    /* Scatter tempv[*] into SPA dense[*] temporarily,
		     * such that tempv[*] can be used for the triangular solve of
		     * the next column of the panel. They will be copied into
		     * ucol[*] after the whole panel has been finished.
		     */
		    isub = lptr + no_zeros;
	/*#pragma ivdep*/
		    for (i = 0; i < segsze; i++) {
			irow = lsub[isub];
			dense_col[dense_col_offset+irow] = tempv[i]; /* Scatter */
			tempv[i] = zero;
			isub++;
	if ( DEBUGlevel>=2 ) {
		if (jj == -1 && krep == 3423)
		    printf("(%d) pdgstrf_bmod1D[scatter] jj %d, dense_col[%d] %e\n",
			   pnum, jj, irow, dense_col[dense_col_offset+irow]);
	}
		    }

		    /* Scatter the update from tempv1[*] into SPA dense[*] */
	/*#pragma ivdep*/
		    for (i = 0; i < nrow; i++) {
			irow = lsub[isub];
	                dense_col[dense_col_offset+irow] -= tempv1[tempv1_offset+i]; /* Scatter-add */
	if (SCATTER_FOUND) {
			if ( col_marker[col_marker_offset+irow] != jj ) {
			    col_marker[col_marker_offset+irow] = jj;
			    col_lsub[col_lsub_offset+w_lsub_end[jj-jcol]++] = irow;
			}
	}
			tempv1[tempv1_offset+i] = zero;
			isub++;
		    }

		} /* else segsze >= 4 ... */

	    } /* for jj ... */

	}

}

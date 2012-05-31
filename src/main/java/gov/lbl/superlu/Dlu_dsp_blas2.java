/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */
/*
 * File name:		dsp_blas2.c
 * Purpose:		Sparse BLAS 2, using some dense BLAS 2 operations.
 */
package gov.lbl.superlu;

import org.netlib.blas.BLAS;

import gov.lbl.superlu.Dlu_supermatrix.NCPformat;
import gov.lbl.superlu.Dlu_supermatrix.NCformat;
import gov.lbl.superlu.Dlu_supermatrix.SCPformat;
import gov.lbl.superlu.Dlu_supermatrix.SuperMatrix;

import static gov.lbl.superlu.Dlu_xerbla_.xerbla_;

import static gov.lbl.superlu.Dlu_pdmemory.doubleCalloc;

import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_ABORT;
import static gov.lbl.superlu.Dlu_slu_mt_util.L_FST_SUPC;
import static gov.lbl.superlu.Dlu_slu_mt_util.L_SUB_START;
import static gov.lbl.superlu.Dlu_slu_mt_util.L_SUB_END;
import static gov.lbl.superlu.Dlu_slu_mt_util.L_LAST_SUPC;
import static gov.lbl.superlu.Dlu_slu_mt_util.L_NZ_START;
import static gov.lbl.superlu.Dlu_slu_mt_util.L_NZ_END;
import static gov.lbl.superlu.Dlu_slu_mt_util.L_SUB;
import static gov.lbl.superlu.Dlu_slu_mt_util.U_NZ_START;
import static gov.lbl.superlu.Dlu_slu_mt_util.U_NZ_END;
import static gov.lbl.superlu.Dlu_slu_mt_util.U_SUB;

import static gov.lbl.superlu.Dlu.USE_VENDOR_BLAS;

import static gov.lbl.superlu.Dlu_dmyblas2.dlsolve;
import static gov.lbl.superlu.Dlu_dmyblas2.dmatvec;
import static gov.lbl.superlu.Dlu_dmyblas2.dusolve;

import static gov.lbl.superlu.Dlu_lsame.lsame_;


public class Dlu_dsp_blas2 {

	static
	int
	sp_dtrsv(char uplo, char trans, char diag, SuperMatrix L,
	         SuperMatrix U, double x[], int[] info)
	{
	/*
	 *   Purpose
	 *   =======
	 *
	 *   sp_dtrsv() solves one of the systems of equations
	 *       A*x = b,   or   A'*x = b,
	 *   where b and x are n element vectors and A is a sparse unit , or
	 *   non-unit, upper or lower triangular matrix.
	 *   No test for singularity or near-singularity is included in this
	 *   routine. Such tests must be performed before calling this routine.
	 *
	 *   Parameters
	 *   ==========
	 *
	 *   uplo   - (input) char*
	 *            On entry, uplo specifies whether the matrix is an upper or
	 *             lower triangular matrix as follows:
	 *                uplo = 'U' or 'u'   A is an upper triangular matrix.
	 *                uplo = 'L' or 'l'   A is a lower triangular matrix.
	 *
	 *   trans  - (input) char*
	 *             On entry, trans specifies the equations to be solved as
	 *             follows:
	 *                trans = 'N' or 'n'   A*x = b.
	 *                trans = 'T' or 't'   A'*x = b.
	 *                trans = 'C' or 'c'   A'*x = b.
	 *
	 *   diag   - (input) char*
	 *             On entry, diag specifies whether or not A is unit
	 *             triangular as follows:
	 *                diag = 'U' or 'u'   A is assumed to be unit triangular.
	 *                diag = 'N' or 'n'   A is not assumed to be unit
	 *                                    triangular.
	 *
	 *   L       - (input) SuperMatrix*
	 *	       The factor L from the factorization Pr*A*Pc=L*U. Use
	 *             compressed row subscripts storage for supernodes,
	 *             i.e., L has types: Stype = SC, Dtype = _D, Mtype = TRLU.
	 *
	 *   U       - (input) SuperMatrix*
	 *	        The factor U from the factorization Pr*A*Pc=L*U.
	 *	        U has types: Stype = NCP, Dtype = _D, Mtype = TRU.
	 *
	 *   x       - (input/output) double*
	 *             Before entry, the incremented array X must contain the n
	 *             element right-hand side vector b. On exit, X is overwritten
	 *             with the solution vector x.
	 *
	 *   info    - (output) int*
	 *             If *info = -i, the i-th argument had an illegal value.
	 *
	 */
	    SCPformat Lstore;
	    NCPformat Ustore;
	    double   Lval[], Uval[];
	    int incx = 1, incy = 1;
	    double alpha = 1.0, beta = 1.0;
	    int fsupc, luptr, istart, irow, k, iptr, jcol, nsuper;
	    int          nsupr, nsupc, nrow, i;
	    double work[];
	    float solve_ops;

	    /* Test the input parameters */
	    info[0] = 0;
	    if ( lsame_(uplo,'L') == 0 && lsame_(uplo, 'U') == 0 ) info[0] = -1;
	    else if ( lsame_(trans, 'N') == 0 && lsame_(trans, 'T') == 0 ) info[0] = -2;
	    else if ( lsame_(diag, 'U') == 0 && lsame_(diag, 'N') == 0 ) info[0] = -3;
	    else if ( L.nrow != L.ncol || L.nrow < 0 ) info[0] = -4;
	    else if ( U.nrow != U.ncol || U.nrow < 0 ) info[0] = -5;
	    if ( info[0] != 0 ) {
		i = -(info[0]);
		xerbla_("sp_dtrsv", i);
		return 0;
	    }

	    Lstore = (SCPformat) L.Store;
	    Lval = (double[]) Lstore.nzval;
	    Ustore = (NCPformat) U.Store;
	    Uval = (double[]) Ustore.nzval;
	    nsuper = Lstore.nsuper;
	    solve_ops = 0;

	    if ( (work = doubleCalloc(L.nrow)) == null )
		SUPERLU_ABORT("Malloc fails for work in sp_dtrsv().");

	    if ( lsame_(trans, 'N') != 0 ) {	/* Form x := inv(A)*x. */

		if ( lsame_(uplo, 'L') != 0 ) {
		    /* Form x := inv(L)*x */
	    	    if ( L.nrow == 0 ) return 0; /* Quick return */

		    for (k = 0; k <= nsuper; k++) {
			fsupc = L_FST_SUPC(Lstore, k);
			istart = L_SUB_START(Lstore, fsupc);
	                nsupr = L_SUB_END(Lstore, fsupc) - istart;
	                nsupc = L_LAST_SUPC(Lstore, k) - fsupc;
			luptr = L_NZ_START(Lstore, fsupc);
			nrow = nsupr - nsupc;

		        solve_ops += nsupc * (nsupc - 1);
		        solve_ops += 2 * nrow * nsupc;

			if ( nsupc == 1 ) {
			    for (iptr=istart+1; iptr < L_SUB_END(Lstore, fsupc); ++iptr) {
				irow = L_SUB(Lstore, iptr);
				++luptr;
				x[irow] -= x[fsupc] * Lval[luptr];
			    }
			} else {
	if (USE_VENDOR_BLAS) {
				BLAS blas = BLAS.getInstance();
			    blas.dtrsv("L", "N", "U", &nsupc, &Lval[luptr], &nsupr,
			       	&x[fsupc], &incx);

			    blas.dgemv("N", &nrow, &nsupc, &alpha, &Lval[luptr+nsupc],
			       	&nsupr, &x[fsupc], &incx, &beta, &work[0], &incy);
	} else {
			    dlsolve (nsupr, nsupc, &Lval[luptr], &x[fsupc]);

			    dmatvec (nsupr, nsupr-nsupc, nsupc, &Lval[luptr+nsupc],
	                             &x[fsupc], &work[0] );
	}

			    iptr = istart + nsupc;
			    for (i = 0; i < nrow; ++i, ++iptr) {
				irow = L_SUB(Lstore, iptr);
				x[irow] -= work[i];	/* Scatter */
				work[i] = 0.0;

			    }
		 	}
		    } /* for k ... */

		} else {
		    /* Form x := inv(U)*x */

		    if ( U.nrow == 0 ) return 0; /* Quick return */

		    for (k = nsuper; k >= 0; k--) {
		    	fsupc = L_FST_SUPC(Lstore, k);
	                nsupr = L_SUB_END(Lstore, fsupc) - L_SUB_START(Lstore, fsupc);
	                nsupc = L_LAST_SUPC(Lstore, k) - fsupc;
		    	luptr = L_NZ_START(Lstore, fsupc);

	    	        solve_ops += nsupc * (nsupc + 1);

			if ( nsupc == 1 ) {
			    x[fsupc] /= Lval[luptr];
			    for (i = U_NZ_START(Ustore, fsupc); i < U_NZ_END(Ustore, fsupc); ++i) {
				irow = U_SUB(Ustore, i);
				x[irow] -= x[fsupc] * Uval[i];
			    }
			} else {
	if (USE_VENDOR_BLAS) {
				BLAS blas = BLAS.getInstance();
			    blas.dtrsv("U", "N", "N", &nsupc, &Lval[luptr], &nsupr,
	                           &x[fsupc], &incx);
	} else {
			    dusolve ( nsupr, nsupc, &Lval[luptr], &x[fsupc] );
	}

	                    for (jcol = fsupc; jcol < fsupc + nsupc; jcol++) {
			        solve_ops += 2*(U_NZ_END(Ustore, jcol) - U_NZ_START(Ustore, jcol));
			    	for (i = U_NZ_START(Ustore, jcol); i < U_NZ_END(Ustore, jcol); i++) {
				    irow = U_SUB(Ustore, i);
				    x[irow] -= x[jcol] * Uval[i];
			    	}
	                    }
			}
		    } /* for k ... */

		}
	    } else { /* Form x := inv(A')*x */

		if ( lsame_(uplo, 'L') != 0 ) {
		    /* Form x := inv(L')*x */
	    	    if ( L.nrow == 0 ) return 0; /* Quick return */

		    for (k = nsuper; k >= 0; --k) {
		    	fsupc = L_FST_SUPC(Lstore, k);
		    	istart = L_SUB_START(Lstore, fsupc);
	                nsupr = L_SUB_END(Lstore, fsupc) - istart;
	                nsupc = L_LAST_SUPC(Lstore, k) - fsupc;
		    	luptr = L_NZ_START(Lstore, fsupc);

			solve_ops += 2 * (nsupr - nsupc) * nsupc;

			for (jcol = fsupc; jcol < L_LAST_SUPC(Lstore, k); jcol++) {
			    iptr = istart + nsupc;
			    for (i = L_NZ_START(Lstore, jcol) + nsupc;
					i < L_NZ_END(Lstore, jcol); i++) {
				irow = L_SUB(Lstore, iptr);
				x[jcol] -= x[irow] * Lval[i];
				iptr++;
			    }
			}

			if ( nsupc > 1 ) {
			    solve_ops += nsupc * (nsupc - 1);

			    BLAS blas = BLAS.getInstance();
			    blas.dtrsv("L", "T", "U", &nsupc, &Lval[luptr], &nsupr,
				&x[fsupc], &incx);
			}
		    }
		} else {
		    /* Form x := inv(U')*x */
		    if ( U.nrow == 0 ) return 0; /* Quick return */

		    for (k = 0; k <= nsuper; k++) {
		    	fsupc = L_FST_SUPC(Lstore, k);
	                nsupr = L_SUB_END(Lstore, fsupc) - L_SUB_START(Lstore, fsupc);
	                nsupc = L_LAST_SUPC(Lstore, k) - fsupc;
		    	luptr = L_NZ_START(Lstore, fsupc);

			for (jcol = fsupc; jcol < fsupc + nsupc; jcol++) {
			    solve_ops += 2*(U_NZ_END(Ustore, jcol) - U_NZ_START(Ustore, jcol));
			    for (i = U_NZ_START(Ustore, jcol); i < U_NZ_END(Ustore, jcol); i++) {
				irow = U_SUB(Ustore, i);
				x[jcol] -= x[irow] * Uval[i];
			    }
			}

			solve_ops += nsupc * (nsupc + 1);

			if ( nsupc == 1 ) {
			    x[fsupc] /= Lval[luptr];
			} else {
				BLAS blas = BLAS.getInstance();
			    blas.dtrsv("U", "T", "N", &nsupc, &Lval[luptr], &nsupr,
				    &x[fsupc], &incx);
			}
		    } /* for k ... */
		}
	    }

	    return 0;
	}



	static
	int
	sp_dgemv(char trans, double alpha, SuperMatrix A, double x[],
		 int incx, double beta, double y[], int incy)
	{
	/*  Purpose
	    =======

	    sp_dgemv()  performs one of the matrix-vector operations
	       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
	    where alpha and beta are scalars, x and y are vectors and A is a
	    sparse A.nrow by A.ncol matrix.

	    Parameters
	    ==========

	    TRANS  - (input) char*
	             On entry, TRANS specifies the operation to be performed as
	             follows:
	                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
	                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
	                TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.

	    ALPHA  - (input) double
	             On entry, ALPHA specifies the scalar alpha.

	    A      - (input) SuperMatrix*
	             Matrix A with a sparse format, of dimension (A.nrow, A.ncol).
	             Currently, the type of A can be:
	                 Stype = NC or NCP; Dtype = SLU_D; Mtype = GE.
	             In the future, more general A can be handled.

	    X      - (input) double*, array of DIMENSION at least
	             ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
	             and at least
	             ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
	             Before entry, the incremented array X must contain the
	             vector x.

	    INCX   - (input) int
	             On entry, INCX specifies the increment for the elements of
	             X. INCX must not be zero.

	    BETA   - (input) double
	             On entry, BETA specifies the scalar beta. When BETA is
	             supplied as zero then Y need not be set on input.

	    Y      - (output) double*,  array of DIMENSION at least
	             ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
	             and at least
	             ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
	             Before entry with BETA non-zero, the incremented array Y
	             must contain the vector y. On exit, Y is overwritten by the
	             updated vector y.

	    INCY   - (input) int
	             On entry, INCY specifies the increment for the elements of
	             Y. INCY must not be zero.

	    ==== Sparse Level 2 Blas routine.
	*/

	    /* Local variables */
	    NCformat Astore;
	    double   Aval[];
	    int info;
	    double temp;
	    int lenx, leny, i, j, irow;
	    int iy, jx, jy, kx, ky;
	    int notran;

	    notran = lsame_(trans, 'N');
	    Astore = (NCformat) A.Store;
	    Aval = Astore.nzval;

	    /* Test the input parameters */
	    info = 0;
	    if ( notran == 0 && lsame_(trans, 'T') == 0 && lsame_(trans, 'C') == 0) info = 1;
	    else if ( A.nrow < 0 || A.ncol < 0 ) info = 3;
	    else if (incx == 0) info = 5;
	    else if (incy == 0)	info = 8;
	    if (info != 0) {
		xerbla_("sp_dgemv ", info);
		return 0;
	    }

	    /* Quick return if possible. */
	    if (A.nrow == 0 || A.ncol == 0 || (alpha == 0. && beta == 1.))
		return 0;

	    /* Set  LENX  and  LENY, the lengths of the vectors x and y, and set
	       up the start points in  X  and  Y. */
	    if (lsame_(trans, 'N') != 0) {
		lenx = A.ncol;
		leny = A.nrow;
	    } else {
		lenx = A.nrow;
		leny = A.ncol;
	    }
	    if (incx > 0) kx = 0;
	    else kx =  - (lenx - 1) * incx;
	    if (incy > 0) ky = 0;
	    else ky =  - (leny - 1) * incy;

	    /* Start the operations. In this version the elements of A are
	       accessed sequentially with one pass through A. */
	    /* First form  y := beta*y. */
	    if (beta != 1.) {
		if (incy == 1) {
		    if (beta == 0.)
			for (i = 0; i < leny; ++i) y[i] = 0.;
		    else
			for (i = 0; i < leny; ++i) y[i] = beta * y[i];
		} else {
		    iy = ky;
		    if (beta == 0.)
			for (i = 0; i < leny; ++i) {
			    y[iy] = 0.;
			    iy += incy;
			}
		    else
			for (i = 0; i < leny; ++i) {
			    y[iy] = beta * y[iy];
			    iy += incy;
			}
		}
	    }

	    if (alpha == 0.) return 0;

	    if ( notran != 0 ) {
		/* Form  y := alpha*A*x + y. */
		jx = kx;
		if (incy == 1) {
		    for (j = 0; j < A.ncol; ++j) {
			if (x[jx] != 0.) {
			    temp = alpha * x[jx];
			    for (i = Astore.colptr[j]; i < Astore.colptr[j+1]; ++i) {
				irow = Astore.rowind[i];
				y[irow] += temp * Aval[i];
			    }
			}
			jx += incx;
		    }
		} else {
		    SUPERLU_ABORT("Not implemented.");
		}
	    } else {
		/* Form  y := alpha*A'*x + y. */
		jy = ky;
		if (incx == 1) {
		    for (j = 0; j < A.ncol; ++j) {
			temp = 0.;
			for (i = Astore.colptr[j]; i < Astore.colptr[j+1]; ++i) {
			    irow = Astore.rowind[i];
			    temp += Aval[i] * x[irow];
			}
			y[jy] += alpha * temp;
			jy += incy;
		    }
		} else {
		    SUPERLU_ABORT("Not implemented.");
		}
	    }
	    return 0;
	} /* sp_dgemv */


}

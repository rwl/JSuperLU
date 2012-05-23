/*! @file Dlu_dgstrs.java
 * \brief Solves a system using LU factorization
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

import lu.jsuper.Dlu_slu_util.SuperLUStat_t;
import lu.jsuper.Dlu_superlu_enum_consts.trans_t;
import lu.jsuper.Dlu_supermatrix.SuperMatrix;

public class Dlu_dgstrs {

	/**! \brief
	 *
	 * <pre>
	 * Purpose
	 * =======
	 *
	 * DGSTRS solves a system of linear equations A*X=B or A'*X=B
	 * with A sparse and B dense, using the LU factorization computed by
	 * DGSTRF.
	 *
	 * See supermatrix.h for the definition of 'SuperMatrix' structure.
	 *
	 * Arguments
	 * =========
	 *
	 * trans   (input) trans_t
	 *          Specifies the form of the system of equations:
	 *          = NOTRANS: A * X = B  (No transpose)
	 *          = TRANS:   A'* X = B  (Transpose)
	 *          = CONJ:    A**H * X = B  (Conjugate transpose)
	 *
	 * L       (input) SuperMatrix*
	 *         The factor L from the factorization Pr*A*Pc=L*U as computed by
	 *         dgstrf(). Use compressed row subscripts storage for supernodes,
	 *         i.e., L has types: Stype = SLU_SC, Dtype = SLU_D, Mtype = SLU_TRLU.
	 *
	 * U       (input) SuperMatrix*
	 *         The factor U from the factorization Pr*A*Pc=L*U as computed by
	 *         dgstrf(). Use column-wise storage scheme, i.e., U has types:
	 *         Stype = SLU_NC, Dtype = SLU_D, Mtype = SLU_TRU.
	 *
	 * perm_c  (input) int*, dimension (L->ncol)
	 *	   Column permutation vector, which defines the
	 *         permutation matrix Pc; perm_c[i] = j means column i of A is
	 *         in position j in A*Pc.
	 *
	 * perm_r  (input) int*, dimension (L->nrow)
	 *         Row permutation vector, which defines the permutation matrix Pr;
	 *         perm_r[i] = j means row i of A is in position j in Pr*A.
	 *
	 * B       (input/output) SuperMatrix*
	 *         B has types: Stype = SLU_DN, Dtype = SLU_D, Mtype = SLU_GE.
	 *         On entry, the right hand side matrix.
	 *         On exit, the solution matrix if info = 0;
	 *
	 * stat     (output) SuperLUStat_t*
	 *          Record the statistics on runtime and floating-point operation count.
	 *          See util.h for the definition of 'SuperLUStat_t'.
	 *
	 * info    (output) int*
	 * 	   = 0: successful exit
	 *	   < 0: if info = -i, the i-th argument had an illegal value
	 * </pre>
	 */
	public static void dgstrs(trans_t trans, SuperMatrix L, SuperMatrix U,
	        int perm_c[], int perm_r[], SuperMatrix B,
	        SuperLUStat_t stat, int info[]) {

	}

}

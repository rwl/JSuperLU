/*! @file Dlu_get_perm_c.java
 * \brief Matrix permutation operations
 *
 * <pre>
 * -- SuperLU routine (version 3.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 1, 2008
 * </pre>
 */
package lu.jsuper;

import lu.jsuper.Dlu_supermatrix.SuperMatrix;


public class Dlu_get_perm_c {

	/**! \brief
	 *
	 * <pre>
	 * Purpose
	 * =======
	 *
	 * GET_PERM_C obtains a permutation matrix Pc, by applying the multiple
	 * minimum degree ordering code by Joseph Liu to matrix A'*A or A+A'.
	 * or using approximate minimum degree column ordering by Davis et. al.
	 * The LU factorization of A*Pc tends to have less fill than the LU
	 * factorization of A.
	 *
	 * Arguments
	 * =========
	 *
	 * ispec   (input) int
	 *         Specifies the type of column ordering to reduce fill:
	 *         = 1: minimum degree on the structure of A^T * A
	 *         = 2: minimum degree on the structure of A^T + A
	 *         = 3: approximate minimum degree for unsymmetric matrices
	 *         If ispec == 0, the natural ordering (i.e., Pc = I) is returned.
	 *
	 * A       (input) SuperMatrix*
	 *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
	 *         of the linear equations is A->nrow. Currently, the type of A
	 *         can be: Stype = NC; Dtype = _D; Mtype = GE. In the future,
	 *         more general A can be handled.
	 *
	 * perm_c  (output) int*
	 *	   Column permutation vector of size A->ncol, which defines the
	 *         permutation matrix Pc; perm_c[i] = j means column i of A is
	 *         in position j in A*Pc.
	 * </pre>
	 */
	public static void get_perm_c(int ispec, SuperMatrix A, int[] perm_c) {

	}

}

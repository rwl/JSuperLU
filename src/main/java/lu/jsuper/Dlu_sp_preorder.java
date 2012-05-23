/*! @file Dlu_sp_preorder.java
 * \brief Permute and performs functions on columns of orginal matrix
 */
package lu.jsuper;

import lu.jsuper.Dlu_slu_util.Dlu_superlu_options_t;
import lu.jsuper.Dlu_supermatrix.SuperMatrix;


public class Dlu_sp_preorder {

	/**! \brief
	 *
	 * <pre>
	 * Purpose
	 * =======
	 *
	 * sp_preorder() permutes the columns of the original matrix. It performs
	 * the following steps:
	 *
	 *    1. Apply column permutation perm_c[] to A's column pointers to form AC;
	 *
	 *    2. If options->Fact = DOFACT, then
	 *       (1) Compute column elimination tree etree[] of AC'AC;
	 *       (2) Post order etree[] to get a postordered elimination tree etree[],
	 *           and a postorder permutation post[];
	 *       (3) Apply post[] permutation to columns of AC;
	 *       (4) Overwrite perm_c[] with the product perm_c * post.
	 *
	 * Arguments
	 * =========
	 *
	 * options (input) superlu_options_t*
	 *         Specifies whether or not the elimination tree will be re-used.
	 *         If options->Fact == DOFACT, this means first time factor A,
	 *         etree is computed, postered, and output.
	 *         Otherwise, re-factor A, etree is input, unchanged on exit.
	 *
	 * A       (input) SuperMatrix*
	 *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
	 *         of the linear equations is A->nrow. Currently, the type of A can be:
	 *         Stype = NC or SLU_NCP; Mtype = SLU_GE.
	 *         In the future, more general A may be handled.
	 *
	 * perm_c  (input/output) int*
	 *	   Column permutation vector of size A->ncol, which defines the
	 *         permutation matrix Pc; perm_c[i] = j means column i of A is
	 *         in position j in A*Pc.
	 *         If options->Fact == DOFACT, perm_c is both input and output.
	 *         On output, it is changed according to a postorder of etree.
	 *         Otherwise, perm_c is input.
	 *
	 * etree   (input/output) int*
	 *         Elimination tree of Pc'*A'*A*Pc, dimension A->ncol.
	 *         If options->Fact == DOFACT, etree is an output argument,
	 *         otherwise it is an input argument.
	 *         Note: etree is a vector of parent pointers for a forest whose
	 *         vertices are the integers 0 to A->ncol-1; etree[root]==A->ncol.
	 *
	 * AC      (output) SuperMatrix*
	 *         The resulting matrix after applied the column permutation
	 *         perm_c[] to matrix A. The type of AC can be:
	 *         Stype = SLU_NCP; Dtype = A->Dtype; Mtype = SLU_GE.
	 * </pre>
	 */
	public static void sp_preorder(Dlu_superlu_options_t options, SuperMatrix A,
			int[] perm_c,  int[] etree, SuperMatrix[] AC) {

	}

}

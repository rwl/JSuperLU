/*! @file Dlu_dgssv.java
 * \brief Solves the system of linear equations A*X=B
 *
 * <pre>
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 * </pre>
 */
package lu.jsuper;

import lu.jsuper.Dlu_slu_util.Dlu_superlu_options_t;
import lu.jsuper.Dlu_slu_util.SuperLUStat_t;
import lu.jsuper.Dlu_superlu_enum_consts.PhaseType;
import lu.jsuper.Dlu_superlu_enum_consts.colperm_t;
import lu.jsuper.Dlu_superlu_enum_consts.fact_t;
import lu.jsuper.Dlu_superlu_enum_consts.trans_t;
import lu.jsuper.Dlu_supermatrix.DNformat;
import lu.jsuper.Dlu_supermatrix.Dtype_t;
import lu.jsuper.Dlu_supermatrix.Mtype_t;
import lu.jsuper.Dlu_supermatrix.Stype_t;
import lu.jsuper.Dlu_supermatrix.SuperMatrix;
import lu.jsuper.Dlu_supermatrix.NRformat;

import static lu.jsuper.Dlu_slu_util.SUPERLU_MAX;
import static lu.jsuper.Dlu_xerbla.xerbla_;
import static lu.jsuper.Dlu_superlu_timer.SuperLU_timer_;
import static lu.jsuper.Dlu_get_perm_c.get_perm_c;
import static lu.jsuper.Dlu_memory.intMalloc;
import static lu.jsuper.Dlu_sp_preorder.sp_preorder;
import static lu.jsuper.Dlu_sp_ienv.sp_ienv;
import static lu.jsuper.Dlu_dutil.dCreate_CompCol_Matrix;
import static lu.jsuper.Dlu_dgstrf.dgstrf;
import static lu.jsuper.Dlu_dgstrs.dgstrs;


public class Dlu_dgssv {

	/**! \brief
	 *
	 * <pre>
	 * Purpose
	 * =======
	 *
	 * DGSSV solves the system of linear equations A*X=B, using the
	 * LU factorization from DGSTRF. It performs the following steps:
	 *
	 *   1. If A is stored column-wise (A.Stype = SLU_NC):
	 *
	 *      1.1. Permute the columns of A, forming A*Pc, where Pc
	 *           is a permutation matrix. For more details of this step,
	 *           see sp_preorder.c.
	 *
	 *      1.2. Factor A as Pr*A*Pc=L*U with the permutation Pr determined
	 *           by Gaussian elimination with partial pivoting.
	 *           L is unit lower triangular with offdiagonal entries
	 *           bounded by 1 in magnitude, and U is upper triangular.
	 *
	 *      1.3. Solve the system of equations A*X=B using the factored
	 *           form of A.
	 *
	 *   2. If A is stored row-wise (A.Stype = SLU_NR), apply the
	 *      above algorithm to the transpose of A:
	 *
	 *      2.1. Permute columns of transpose(A) (rows of A),
	 *           forming transpose(A)*Pc, where Pc is a permutation matrix.
	 *           For more details of this step, see sp_preorder.c.
	 *
	 *      2.2. Factor A as Pr*transpose(A)*Pc=L*U with the permutation Pr
	 *           determined by Gaussian elimination with partial pivoting.
	 *           L is unit lower triangular with offdiagonal entries
	 *           bounded by 1 in magnitude, and U is upper triangular.
	 *
	 *      2.3. Solve the system of equations A*X=B using the factored
	 *           form of A.
	 *
	 *   See supermatrix.h for the definition of 'SuperMatrix' structure.
	 *
	 * Arguments
	 * =========
	 *
	 * options (input) superlu_options_t*
	 *         The structure defines the input parameters to control
	 *         how the LU decomposition will be performed and how the
	 *         system will be solved.
	 *
	 * A       (input) SuperMatrix*
	 *         Matrix A in A*X=B, of dimension (A.nrow, A.ncol). The number
	 *         of linear equations is A.nrow. Currently, the type of A can be:
	 *         Stype = SLU_NC or SLU_NR; Dtype = SLU_D; Mtype = SLU_GE.
	 *         In the future, more general A may be handled.
	 *
	 * perm_c  (input/output) int*
	 *         If A.Stype = SLU_NC, column permutation vector of size A.ncol
	 *         which defines the permutation matrix Pc; perm_c[i] = j means
	 *         column i of A is in position j in A*Pc.
	 *         If A.Stype = SLU_NR, column permutation vector of size A.nrow
	 *         which describes permutation of columns of transpose(A)
	 *         (rows of A) as described above.
	 *
	 *         If options.ColPerm = MY_PERMC or options.Fact = SamePattern or
	 *            options.Fact = SamePattern_SameRowPerm, it is an input argument.
	 *            On exit, perm_c may be overwritten by the product of the input
	 *            perm_c and a permutation that postorders the elimination tree
	 *            of Pc'*A'*A*Pc; perm_c is not changed if the elimination tree
	 *            is already in postorder.
	 *         Otherwise, it is an output argument.
	 *
	 * perm_r  (input/output) int*
	 *         If A.Stype = SLU_NC, row permutation vector of size A.nrow,
	 *         which defines the permutation matrix Pr, and is determined
	 *         by partial pivoting.  perm_r[i] = j means row i of A is in
	 *         position j in Pr*A.
	 *         If A.Stype = SLU_NR, permutation vector of size A.ncol, which
	 *         determines permutation of rows of transpose(A)
	 *         (columns of A) as described above.
	 *
	 *         If options.RowPerm = MY_PERMR or
	 *            options.Fact = SamePattern_SameRowPerm, perm_r is an
	 *            input argument.
	 *         otherwise it is an output argument.
	 *
	 * L       (output) SuperMatrix*
	 *         The factor L from the factorization
	 *             Pr*A*Pc=L*U              (if A.Stype = SLU_NC) or
	 *             Pr*transpose(A)*Pc=L*U   (if A.Stype = SLU_NR).
	 *         Uses compressed row subscripts storage for supernodes, i.e.,
	 *         L has types: Stype = SLU_SC, Dtype = SLU_D, Mtype = SLU_TRLU.
	 *
	 * U       (output) SuperMatrix*
	 *	   The factor U from the factorization
	 *             Pr*A*Pc=L*U              (if A.Stype = SLU_NC) or
	 *             Pr*transpose(A)*Pc=L*U   (if A.Stype = SLU_NR).
	 *         Uses column-wise storage scheme, i.e., U has types:
	 *         Stype = SLU_NC, Dtype = SLU_D, Mtype = SLU_TRU.
	 *
	 * B       (input/output) SuperMatrix*
	 *         B has types: Stype = SLU_DN, Dtype = SLU_D, Mtype = SLU_GE.
	 *         On entry, the right hand side matrix.
	 *         On exit, the solution matrix if info = 0;
	 *
	 * stat   (output) SuperLUStat_t*
	 *        Record the statistics on runtime and floating-point operation count.
	 *        See util.h for the definition of 'SuperLUStat_t'.
	 *
	 * info    (output) int*
	 *	   = 0: successful exit
	 *         > 0: if info = i, and i is
	 *             <= A.ncol: U(i,i) is exactly zero. The factorization has
	 *                been completed, but the factor U is exactly singular,
	 *                so the solution could not be computed.
	 *             > A.ncol: number of bytes allocated when memory allocation
	 *                failure occurred, plus A.ncol.
	 * </pre>
	 */
	public static void dgssv(Dlu_superlu_options_t options, SuperMatrix A,
			int perm_c[], int perm_r[],
			SuperMatrix[] L, SuperMatrix[] U, SuperMatrix B,
			SuperLUStat_t stat, int[] info) {

	    DNformat Bstore;
	    /* A in SLU_NC format used by the factorization routine.*/
	    SuperMatrix AA = null;
	    /* Matrix postmultiplied by Pc */
	    SuperMatrix[] AC = new SuperMatrix[1];
	    int      lwork = 0, etree[], i;

	    /* Set default values for some parameters */
	    int      panel_size;     /* panel size */
	    int      relax;          /* no of columns in a relaxed snodes */
	    int      permc_spec;
	    trans_t  trans = trans_t.NOTRANS;
	    double   utime[];
	    double   t;	/* Temporary time */

	    /* Test the input parameters ... */
	    info[0] = 0;
	    Bstore = (DNformat) B.Store;
	    if ( options.Fact != fact_t.DOFACT ) info[0] = -1;
	    else if ( A.nrow != A.ncol || A.nrow < 0 ||
		 (A.Stype != Stype_t.SLU_NC && A.Stype != Stype_t.SLU_NR) ||
		 A.Dtype != Dtype_t.SLU_D || A.Mtype != Mtype_t.SLU_GE )
		info[0] = -2;
	    else if ( B.ncol < 0 || Bstore.lda < SUPERLU_MAX(0, A.nrow) ||
		B.Stype != Stype_t.SLU_DN || B.Dtype != Dtype_t.SLU_D ||
		B.Mtype != Mtype_t.SLU_GE )
		info[0] = -7;
	    if ( info[0] != 0 ) {
		i = -(info[0]);
		xerbla_("dgssv", i);
		return;
	    }

	    utime = stat.utime;

	    /* Convert A to SLU_NC format when necessary. */
	    if ( A.Stype == Stype_t.SLU_NR ) {
		NRformat Astore = (NRformat) A.Store;
		AA = dCreate_CompCol_Matrix(A.ncol, A.nrow, Astore.nnz,
				       Astore.nzval, Astore.colind, Astore.rowptr,
				       Stype_t.SLU_NC, A.Dtype, A.Mtype);
		trans = trans_t.TRANS;
	    } else {
	        if ( A.Stype == Stype_t.SLU_NC ) AA = A;
	    }

	    t = SuperLU_timer_();
	    /*
	     * Get column permutation vector perm_c[], according to permc_spec:
	     *   permc_spec = NATURAL:  natural ordering
	     *   permc_spec = MMD_AT_PLUS_A: minimum degree on structure of A'+A
	     *   permc_spec = MMD_ATA:  minimum degree on structure of A'*A
	     *   permc_spec = COLAMD:   approximate minimum degree column ordering
	     *   permc_spec = MY_PERMC: the ordering already supplied in perm_c[]
	     */
	    permc_spec = options.ColPerm.ordinal();
	    if ( permc_spec != colperm_t.MY_PERMC.ordinal() &&
	    		options.Fact == fact_t.DOFACT )
	      get_perm_c(permc_spec, AA, perm_c);
	    utime[PhaseType.COLPERM.ordinal()] = SuperLU_timer_() - t;

	    etree = intMalloc(A.ncol);

	    t = SuperLU_timer_();
	    sp_preorder(options, AA, perm_c, etree, AC);
	    utime[PhaseType.ETREE.ordinal()] = SuperLU_timer_() - t;

	    panel_size = sp_ienv(1);
	    relax = sp_ienv(2);

	    /*printf("Factor PA = LU ... relax %d\tw %d\tmaxsuper %d\trowblk %d\n",
		  relax, panel_size, sp_ienv(3), sp_ienv(4));*/
	    t = SuperLU_timer_();
	    /* Compute the LU factorization of A. */
	    dgstrf(options, AC[0], relax, panel_size, etree,
	            null, lwork, perm_c, perm_r, L, U, stat, info);
	    utime[PhaseType.FACT.ordinal()] = SuperLU_timer_() - t;

	    t = SuperLU_timer_();
	    if ( info[0] == 0 ) {
	        /* Solve the system A*X=B, overwriting B with X. */
	        dgstrs (trans, L[0], U[0], perm_c, perm_r, B, stat, info);
	    }
	    utime[PhaseType.SOLVE.ordinal()] = SuperLU_timer_() - t;

	    etree = null;
	    AC = null;
	    if ( A.Stype == Stype_t.SLU_NR ) {
		AA = null;
	    }

	}

}

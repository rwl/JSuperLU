/** @file Dlu_slu_util.java
 * \brief Utility header file
 *
 * -- SuperLU routine (version 4.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November, 2010
 *
 */
package lu.jsuper;

import lu.jsuper.Dlu_superlu_enum_consts.IterRefine_t;
import lu.jsuper.Dlu_superlu_enum_consts.colperm_t;
import lu.jsuper.Dlu_superlu_enum_consts.fact_t;
import lu.jsuper.Dlu_superlu_enum_consts.milu_t;
import lu.jsuper.Dlu_superlu_enum_consts.norm_t;
import lu.jsuper.Dlu_superlu_enum_consts.rowperm_t;
import lu.jsuper.Dlu_superlu_enum_consts.trans_t;
import lu.jsuper.Dlu_superlu_enum_consts.yes_no_t;

import static lu.jsuper.Dlu_util.superlu_abort_and_exit;


public class Dlu_slu_util {

	public static void USER_ABORT(String msg) {
		superlu_abort_and_exit(msg);
	}

	public static void ABORT(String err_msg) {
		ABORT(err_msg, null);
	}

	public static void ABORT(String err_msg, Throwable e) {
		String msg;
		if (e == null) {
			msg = err_msg;
		} else {
			msg = String.format("%s at line %d in class %s\n", err_msg,
					e.getStackTrace()[2].getLineNumber(),
					e.getStackTrace()[2].getClassName());
		}
		USER_ABORT(msg);
	}

	public static int SUPERLU_MAX(int x, int y) {
		return (x) > (y) ? (x) : (y);
	}

	public static int SUPERLU_MIN(int x, int y) {
		return (x) < (y) ? (x) : (y);
	}

	/**
	 *-- This contains the options used to control the solution process.
	 *
	 * Fact   (fact_t)
	 *        Specifies whether or not the factored form of the matrix
	 *        A is supplied on entry, and if not, how the matrix A should
	 *        be factorizaed.
	 *        = DOFACT: The matrix A will be factorized from scratch, and the
	 *             factors will be stored in L and U.
	 *        = SamePattern: The matrix A will be factorized assuming
	 *             that a factorization of a matrix with the same sparsity
	 *             pattern was performed prior to this one. Therefore, this
	 *             factorization will reuse column permutation vector
	 *             ScalePermstruct->perm_c and the column elimination tree
	 *             LUstruct->etree.
	 *        = SamePattern_SameRowPerm: The matrix A will be factorized
	 *             assuming that a factorization of a matrix with the same
	 *             sparsity	pattern and similar numerical values was performed
	 *             prior to this one. Therefore, this factorization will reuse
	 *             both row and column scaling factors R and C, both row and
	 *             column permutation vectors perm_r and perm_c, and the
	 *             data structure set up from the previous symbolic factorization.
	 *        = FACTORED: On entry, L, U, perm_r and perm_c contain the
	 *              factored form of A. If DiagScale is not NOEQUIL, the matrix
	 *              A has been equilibrated with scaling factors R and C.
	 *
	 * Equil  (yes_no_t)
	 *        Specifies whether to equilibrate the system (scale A's row and
	 *        columns to have unit norm).
	 *
	 * ColPerm (colperm_t)
	 *        Specifies what type of column permutation to use to reduce fill.
	 *        = NATURAL: use the natural ordering
	 *        = MMD_ATA: use minimum degree ordering on structure of A'*A
	 *        = MMD_AT_PLUS_A: use minimum degree ordering on structure of A'+A
	 *        = COLAMD: use approximate minimum degree column ordering
	 *        = MY_PERMC: use the ordering specified by the user
	 *
	 * Trans  (trans_t)
	 *        Specifies the form of the system of equations:
	 *        = NOTRANS: A * X = B        (No transpose)
	 *        = TRANS:   A**T * X = B     (Transpose)
	 *        = CONJ:    A**H * X = B     (Transpose)
	 *
	 * IterRefine (IterRefine_t)
	 *        Specifies whether to perform iterative refinement.
	 *        = NO: no iterative refinement
	 *        = SLU_SINGLE: perform iterative refinement in single precision
	 *        = SLU_DOUBLE: perform iterative refinement in double precision
	 *        = SLU_EXTRA: perform iterative refinement in extra precision
	 *
	 * DiagPivotThresh (double, in [0.0, 1.0]) (only for sequential SuperLU)
	 *        Specifies the threshold used for a diagonal entry to be an
	 *        acceptable pivot.
	 *
	 * SymmetricMode (yest_no_t)
	 *        Specifies whether to use symmetric mode. Symmetric mode gives
	 *        preference to diagonal pivots, and uses an (A'+A)-based column
	 *        permutation algorithm.
	 *
	 * PivotGrowth (yes_no_t)
	 *        Specifies whether to compute the reciprocal pivot growth.
	 *
	 * ConditionNumber (ues_no_t)
	 *        Specifies whether to compute the reciprocal condition number.
	 *
	 * RowPerm (rowperm_t) (only for SuperLU_DIST or ILU)
	 *        Specifies whether to permute rows of the original matrix.
	 *        = NO: not to permute the rows
	 *        = LargeDiag: make the diagonal large relative to the off-diagonal
	 *        = MY_PERMR: use the permutation given by the user
	 *
	 * ILU_DropRule (int)
	 *        Specifies the dropping rule:
	 *	  = DROP_BASIC:   Basic dropping rule, supernodal based ILUTP(tau).
	 *	  = DROP_PROWS:   Supernodal based ILUTP(p,tau), p = gamma * nnz(A)/n.
	 *	  = DROP_COLUMN:  Variant of ILUTP(p,tau), for j-th column,
	 *			      p = gamma * nnz(A(:,j)).
	 *	  = DROP_AREA:    Variation of ILUTP, for j-th column, use
	 *			      nnz(F(:,1:j)) / nnz(A(:,1:j)) to control memory.
	 *	  = DROP_DYNAMIC: Modify the threshold tau during factorizaion:
	 *			  If nnz(L(:,1:j)) / nnz(A(:,1:j)) > gamma
	 *				  tau_L(j) := MIN(tau_0, tau_L(j-1) * 2);
	 *			  Otherwise
	 *				  tau_L(j) := MAX(tau_0, tau_L(j-1) / 2);
	 *			  tau_U(j) uses the similar rule.
	 *			  NOTE: the thresholds used by L and U are separate.
	 *	  = DROP_INTERP:  Compute the second dropping threshold by
	 *	                  interpolation instead of sorting (default).
	 *  		          In this case, the actual fill ratio is not
	 *			  guaranteed to be smaller than gamma.
	 *   	  Note: DROP_PROWS, DROP_COLUMN and DROP_AREA are mutually exclusive.
	 *	  ( Default: DROP_BASIC | DROP_AREA )
	 *
	 * ILU_DropTol (double)
	 *        numerical threshold for dropping.
	 *
	 * ILU_FillFactor (double)
	 *        Gamma in the secondary dropping.
	 *
	 * ILU_Norm (norm_t)
	 *        Specify which norm to use to measure the row size in a
	 *        supernode: infinity-norm, 1-norm, or 2-norm.
	 *
	 * ILU_FillTol (double)
	 *        numerical threshold for zero pivot perturbation.
	 *
	 * ILU_MILU (milu_t)
	 *        Specifies which version of MILU to use.
	 *
	 * ILU_MILU_Dim (double)
	 *        Dimension of the PDE if available.
	 *
	 * ReplaceTinyPivot (yes_no_t) (only for SuperLU_DIST)
	 *        Specifies whether to replace the tiny diagonals by
	 *        sqrt(epsilon)*||A|| during LU factorization.
	 *
	 * SolveInitialized (yes_no_t) (only for SuperLU_DIST)
	 *        Specifies whether the initialization has been performed to the
	 *        triangular solve.
	 *
	 * RefineInitialized (yes_no_t) (only for SuperLU_DIST)
	 *        Specifies whether the initialization has been performed to the
	 *        sparse matrix-vector multiplication routine needed in iterative
	 *        refinement.
	 *
	 * PrintStat (yes_no_t)
	 *        Specifies whether to print the solver's statistics.
	 */
	public static class Dlu_superlu_options_t {
	    public fact_t        Fact;
	    public yes_no_t      Equil;
	    public colperm_t     ColPerm;
	    public trans_t       Trans;
	    public IterRefine_t  IterRefine;
	    public double        DiagPivotThresh;
	    public yes_no_t      SymmetricMode;
	    public yes_no_t      PivotGrowth;
	    public yes_no_t      ConditionNumber;
	    public rowperm_t     RowPerm;
	    public int 	  ILU_DropRule;
	    public double	  ILU_DropTol;    /* threshold for dropping */
	    public double	  ILU_FillFactor; /* gamma in the secondary dropping */
	    public norm_t	  ILU_Norm;       /* infinity-norm, 1-norm, or 2-norm */
	    public double	  ILU_FillTol;    /* threshold for zero pivot perturbation */
	    public milu_t	  ILU_MILU;
	    public double	  ILU_MILU_Dim;   /* Dimension of PDE (if available) */
	    public yes_no_t      ParSymbFact;
	    public yes_no_t      ReplaceTinyPivot; /* used in SuperLU_DIST */
	    public yes_no_t      SolveInitialized;
	    public yes_no_t      RefineInitialized;
	    public yes_no_t      PrintStat;
	    public int           nnzL, nnzU;      /* used to store nnzs for now       */
	    public int           num_lookaheads;  /* num of levels in look-ahead      */
	    public yes_no_t      lookahead_etree; /* use etree computed from the
					      serial symbolic factorization */
	    public yes_no_t      SymPattern;      /* symmetric factorization          */
	}

	public static class SuperLUStat_t {
	    public int     panel_histo[]; /* histogram of panel size distribution */
	    public double  utime[];       /* running time at various phases */
	    public float   ops[];         /* operation count at various phases */
	    public int     TinyPivots;    /* number of tiny pivots */
	    public int     RefineSteps;   /* number of iterative refinement steps */
	    public int     expansions;    /* number of memory expansions */
	}
}

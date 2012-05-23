/* @file Dlu_superlu_enum_consts.java
 * \brief enum constants header file
 *
 * -- SuperLU routine (version 4.1) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * October 1, 2010
 *
 */
package lu.jsuper;


public class Dlu_superlu_enum_consts {

	/***********************************************************************
	 * Enumerate types
	 ***********************************************************************/
	public enum yes_no_t {NO, YES}
	public enum fact_t {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED}
	public enum rowperm_t {NOROWPERM, LargeDiag, MY_PERMR}
	public enum colperm_t {NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD,
	      METIS_AT_PLUS_A, PARMETIS, ZOLTAN, MY_PERMC}
	public enum trans_t {NOTRANS, TRANS, CONJ}
	public enum DiagScale_t {NOEQUIL, ROW, COL, BOTH}
	public enum IterRefine_t {NOREFINE, SLU_SINGLE_1, SLU_DOUBLE, SLU_EXTRA}
	public enum MemType {LUSUP, UCOL, LSUB, USUB, LLVL, ULVL}
	public enum stack_end_t {HEAD, TAIL}
	public enum LU_space_t {SYSTEM, USER}
	public enum norm_t {ONE_NORM, TWO_NORM, INF_NORM}
	public enum milu_t {SILU, SMILU_1, SMILU_2, SMILU_3}

	/**
	 * The following enumerate type is used by the statistics variable
	 * to keep track of flop count and time spent at various stages.
	 *
	 * Note that not all of the fields are disjoint.
	 */
	public enum PhaseType {
	    COLPERM, /* find a column ordering that minimizes fills */
	    ROWPERM, /* find a row ordering maximizes diagonal. */
	    RELAX,   /* find artificial supernodes */
	    ETREE,   /* compute column etree */
	    EQUIL,   /* equilibrate the original matrix */
	    SYMBFAC, /* symbolic factorization. */
	    DIST,    /* distribute matrix. */
	    FACT,    /* perform LU factorization */
	    COMM,    /* communication for factorization */
	    SOL_COMM,/* communication for solve */
	    RCOND,   /* estimate reciprocal condition number */
	    SOLVE,   /* forward and back solves */
	    REFINE,  /* perform iterative refinement */
	    TRSV,    /* fraction of FACT spent in xTRSV */
	    GEMV,    /* fraction of FACT spent in xGEMV */
	    FERR,    /* estimate error bounds after iterative refinement */
	    NPHASES  /* total number of phases */
	}

}

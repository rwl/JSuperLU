/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 */
package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.pxgstrf_shared_t;
import gov.lbl.superlu.Dlu_supermatrix.NCPformat;
import gov.lbl.superlu.Dlu_supermatrix.SCPformat;

import static gov.lbl.superlu.Dlu.PROFILE;
import static gov.lbl.superlu.Dlu_superlu_timer.SuperLU_timer_;
import static gov.lbl.superlu.Dlu_util.superlu_abort_and_exit;



public class Dlu_slu_mt_util {

	public static void USER_ABORT(String msg) {
		superlu_abort_and_exit(msg);
	}

	public static void SUPERLU_ABORT(String err_msg) {
		SUPERLU_ABORT(err_msg, null);
	}

	public static void SUPERLU_ABORT(String err_msg, Throwable e) {
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

	public static int MAX(int x, int y) {
		return (x) > (y) ? (x) : (y);
	}

	public static int MIN(int x, int y) {
		return (x) < (y) ? (x) : (y);
	}

	public static double MAX(double x, double y) {
		return (x) > (y) ? (x) : (y);
	}

	public static double MIN(double x, double y) {
		return (x) < (y) ? (x) : (y);
	}

	public static int SUPERLU_MAX(int x, int y) {
		return (x) > (y) ? (x) : (y);
	}

	public static int SUPERLU_MIN(int x, int y) {
		return (x) < (y) ? (x) : (y);
	}

	public static double SUPERLU_MAX(double x, double y) {
		return (x) > (y) ? (x) : (y);
	}

	public static double SUPERLU_MIN(double x, double y) {
		return (x) < (y) ? (x) : (y);
	}

	/*********************************************************
	 * Macros used for easy access of sparse matrix entries. *
	 *********************************************************/
	static int L_SUB_START(SCPformat Lstore, int col) {
		return Lstore.rowind_colbeg[col];
	}
	static int L_SUB_END(SCPformat Lstore, int col) {
		return Lstore.rowind_colend[col];
	}
	static int L_SUB(SCPformat Lstore, int ptr) {
		return Lstore.rowind[ptr];
	}
	static int L_NZ_START(SCPformat Lstore, int col) {
		return Lstore.nzval_colbeg[col];
	}
	static int L_NZ_END(SCPformat Lstore, int col) {
		return Lstore.nzval_colend[col];
	}
	static int L_FST_SUPC(SCPformat Lstore, int superno) {
		return Lstore.sup_to_colbeg[superno];
	}
	static int L_LAST_SUPC(SCPformat Lstore, int superno) {
		return Lstore.sup_to_colend[superno];
	}
	static int U_NZ_START(NCPformat Ustore, int col) {
		return Ustore.colbeg[col];
	}
	static int U_NZ_END(NCPformat Ustore, int col) {
		return Ustore.colend[col];
	}
	static int U_SUB(NCPformat Ustore, int ptr) {
		return Ustore.rowind[ptr];
	}

	static int SUPER_REP(int[] xsup_end, int s) {
		return xsup_end[s]-1;
	}
	static int SUPER_FSUPC(int[] xsup, int s) {
		return xsup[s];
	}
	static boolean SINGLETON(int[] xsup_end, int[] xsup, int s) {
		return (xsup_end[s] - xsup[s]) == 1;
	}
	static int ISPRUNED(int[] ispruned, int j) {
		return ispruned[j];
	}
	static int STATE(pxgstrf_shared_t pxgstrf_shared, int j) {
		return pxgstrf_shared.pan_status[j].state.ordinal();
	}
	static int DADPANEL(int[] etree, pxgstrf_shared_t pxgstrf_shared, int j) {
		return etree[j + pxgstrf_shared.pan_status[j].size-1];
	}

	static void TIC(double[] t) {
		if (PROFILE) {
			t[0] = SuperLU_timer_();
		}
	}

	static void TOC(double t2[], double t1) {
		if (PROFILE) {
			t2[0] = SuperLU_timer_() - t1;
		}
	}

	/*
	 * Constants
	 */
	public static final int EMPTY =	(-1);
	public static final int FALSE =	0;
	public static final int TRUE =	1;

	/**********************
	  Enumerated constants
	  *********************/
	public enum yes_no_t {NO, YES}
	public enum trans_t {NOTRANS, TRANS, CONJ}
	public enum fact_t {DOFACT, EQUILIBRATE, FACTORED}
	public enum colperm_t {NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD,
		      METIS_AT_PLUS_A, PARMETIS, MY_PERMC}
	public enum equed_t {NOEQUIL, ROW, COL, BOTH}
	public enum MemType {LUSUP, UCOL, LSUB, USUB}

	/* Number of marker arrays used in the symbolic factorization,
	   each of size nrow. */
	public static final int NO_MARKER =     3;

	public static final int LOCOL =    70;
	public static final int HICOL =    78;
	public static final int BADROW =   44;
	public static final int BADCOL =   35;
	public static final int BADPAN =   BADCOL;
	public static final int BADREP =   35;

	public enum PhaseType {
	    RELAX,
	    COLPERM,
	    ETREE,
	    EQUIL,
	    FINDDOMAIN,
	    FACT,
	    DFS,
	    FLOAT,
	    TRSV,
	    GEMV,
	    RCOND,
	    TRISOLVE,
	    SOLVE,
	    REFINE,
	    FERR,
	    NPHASES
	}

	/*
	 * *********************************************************************
	 * The superlumt_options_t structure contains the shared variables used
	 * for factorization, which are passed to each thread.
	 * *********************************************************************
	 *
	 * nprocs (int)
	 *        Number of processes (or threads) to be spawned and used to perform
	 *        the LU factorization by pdgstrf().
	 *
	 * fact   (fact_t)
	 *        Specifies whether or not the factored form of the matrix
	 *        A is supplied on entry, and if not, whether the matrix A should
	 *        be equilibrated before it is factored.
	 *        = DOFACT: The matrix A will be factored, and the factors will be
	 *              stored in L and U.
	 *        = EQUILIBRATE: The matrix A will be equilibrated if necessary, then
	 *              factored into L and U.
	 *
	 * trans  (trans_t)
	 *        Specifies the form of the system of equations:
	 *        = NOTRANS: A * X = B        (No transpose)
	 *        = TRANS:   A**T * X = B     (Transpose)
	 *        = CONJ:    A**H * X = B     (Transpose)
	 *
	 * refact (yes_no_t)
	 *        Specifies whether this is first time or subsequent factorization.
	 *        = NO:  this factorization is treated as the first one;
	 *        = YES: it means that a factorization was performed prior to this
	 *               one. Therefore, this factorization will re-use some
	 *               existing data structures, such as L and U storage, column
	 *               elimination tree, and the symbolic information of the
	 *               Householder matrix.
	 *
	 * panel_size (int)
	 *        A panel consists of at most panel_size consecutive columns.
	 *
	 * relax  (int)
	 *        To control degree of relaxing supernodes. If the number
	 *        of nodes (columns) in a subtree of the elimination tree is less
	 *        than relax, this subtree is considered as one supernode,
	 *        regardless of the row structures of those columns.
	 *
	 * diag_pivot_thresh (double)
	 *        Diagonal pivoting threshold. At step j of the Gaussian elimination,
	 *        if abs(A_jj) >= diag_pivot_thresh * (max_(i>=j) abs(A_ij)),
	 *        use A_jj as pivot, else use A_ij with maximum magnitude.
	 *        0 <= diag_pivot_thresh <= 1. The default value is 1,
	 *        corresponding to partial pivoting.
	 *
	 * drop_tol (double) (NOT IMPLEMENTED)
	 *	  Drop tolerance parameter. At step j of the Gaussian elimination,
	 *        if abs(A_ij)/(max_i abs(A_ij)) < drop_tol, drop entry A_ij.
	 *        0 <= drop_tol <= 1. The default value of drop_tol is 0,
	 *        corresponding to not dropping any entry.
	 *
	 * usepr  (yes_no_t)
	 *        Whether the pivoting will use perm_r specified by the user.
	 *        = YES: use perm_r; perm_r is input, unchanged on exit.
	 *        = NO:  perm_r is determined by partial pivoting, and is output.
	 *
	 * SymmetricMode (yest_no_t)
	 *        Specifies whether to use symmetric mode.
	 *
	 * PrintStat (yes_no_t)
	 *        Specifies whether to print solver's statistics.
	 *
	 * perm_c (int*) dimension A->ncol
	 *	  Column permutation vector, which defines the
	 *        permutation matrix Pc; perm_c[i] = j means column i of A is
	 *        in position j in A*Pc.
	 *        When search for diagonal, perm_c[*] is applied to the
	 *        row subscripts of A, so that diagonal threshold pivoting
	 *        can find the diagonal of A, instead of that of A*Pc.
	 *
	 * perm_r (int*) dimension A->nrow
	 *        Row permutation vector which defines the permutation matrix Pr,
	 *        perm_r[i] = j means row i of A is in position j in Pr*A.
	 *        If usepr = NO, perm_r is output argument;
	 *        If usepr = YES, the pivoting routine will try to use the input
	 *           perm_r, unless a certain threshold criterion is violated.
	 *           In that case, perm_r is overwritten by a new permutation
	 *           determined by partial pivoting or diagonal threshold pivoting.
	 *
	 * work   (void*) of size lwork
	 *        User-supplied work space and space for the output data structures.
	 *        Not referenced if lwork = 0;
	 *
	 * lwork  (int)
	 *        Specifies the length of work array.
	 *        = 0:  allocate space internally by system malloc;
	 *        > 0:  use user-supplied work array of length lwork in bytes,
	 *              returns error if space runs out.
	 *        = -1: the routine guesses the amount of space needed without
	 *              performing the factorization, and returns it in
	 *              superlu_memusage->total_needed; no other side effects.
	 *
	 * etree  (int*)
	 *        Elimination tree of A'*A, dimension A->ncol.
	 *        Note: etree is a vector of parent pointers for a forest whose
	 *        vertices are the integers 0 to A->ncol-1; etree[root]==A->ncol.
	 *        On input, the columns of A should be permutated so that the
	 *        etree is in a certain postorder.
	 *
	 * colcnt_h (int*)
	 *        Column colunts of the Householder matrix.
	 *
	 * part_super_h (int*)
	 *        Partition of the supernodes in the Householder matrix.
	 *	  part_super_h[k] = size of the supernode beginning at column k;
	 * 	                  = 0, elsewhere.
	 *
	 *
	 */
	public static class superlumt_options_t {
		public int        nprocs;
		public fact_t     fact;
		public trans_t    trans;
		public yes_no_t   refact;
		public int        panel_size;
		public int        relax;
		public double     diag_pivot_thresh;
		public double     drop_tol;
		public colperm_t  ColPerm;
		public yes_no_t[]   usepr = new yes_no_t[1];
		public yes_no_t   SymmetricMode;
		public yes_no_t   PrintStat;

	    /* The following arrays are persistent during repeated factorizations. */
		public int  perm_c[];
		public int  perm_r[];
		public double work[];
		public int  lwork;

	    /* The following structural arrays are computed internally by
	       dsp_colorder(), so the user does not provide them on input.
	       These 3 arrays are computed in the first factorization, and are
	       re-used in the subsequent factors of the matrices with the same
	       nonzero structure. */
		public int  etree[];
		public int  colcnt_h[];
		public int  part_super_h[];
	}

	/* ----------------------------------------------
	    The definitions below are used for profiling.
	   ---------------------------------------------- */

	/* The statistics to be kept by each processor. */
	public static class procstat_t {
		public int	    panels;    /* number of panels taken */
		public float   fcops;     /* factor floating-point operations */
		public double  fctime;    /* factor time */
		public int     skedwaits; /* how many times the processor fails to get a task */
		public double  skedtime;  /* time spent in the scheduler */
		public double  cs_time;   /* time spent in the critical sections */
		public double  spintime;  /* spin-wait time */
		public int     pruned;
		public int     unpruned;
	}


	/* Statistics about each panel. */

	public static class panstat_t {
		public int    size;      /* size of the panel */
		public int    pnum;      /* which processor grabs this panel */
		public double starttime; /* at waht time this panel is assigned to a proc */
		public double fctime;    /* factorization time */
		public float  flopcnt;   /* floating-point operations */
		public int    pipewaits; /* how many times the panel waited during pipelining */
		public double spintime;  /* spin waiting time during pipelining */
	}

	/* How was a panel selected by the scheduler */
	public enum how_selected_t {NOPIPE, DADPAN, PIPE}

	/* Headers for 4 types of dynamatically managed memory */
	public static class ExpHeader {//e_node {
		public int size;      /* length of the memory that has been used */
		public double mem[];  /* pointer to the new malloc'd store */
	}

	/* The structure to keep track of memory usage. */
	public static class superlu_memusage_t {
		public float for_lu;
		public float total_needed;
		public int   expansions;
	}

	public static class stat_relax_t {
		public float flops;
		public int     nzs;
		public double  fctime;
	}

	public static class stat_col_t {
		public float flops;
		public int nzs;
		public double fctime;
	}

	public static class stat_snode_t {
		public int ncols;
		public float flops;
		public int nzs;
		public double fctime;
	}

	/* -------------------------------------------------------------------
	   The definitions below are used to simulate parallel execution time.
	   ------------------------------------------------------------------- */
	public static class cp_panel_t {
		public float est;  /* earliest (possible) start time of the panel */
		public float pdiv; /* time in flops spent in the (inner) panel factorization */
	}

	public static class desc_eft_t {
		public float eft;  /* earliest finishing time */
		public float pmod; /* pmod update to the ancestor panel */
	}

	/* All statistics. */
	public static class Gstat_t {
		public int     	panel_histo[];	/* Panel size distribution */
		public double  	utime[];
		public float 	ops[];
		public procstat_t 	procstat[];
		public panstat_t	panstat[];
		public int      	num_panels;
		public float     	dom_flopcnt;
		public float     	flops_last_P_panels;
	    /**/
		public stat_relax_t stat_relax[];
		public stat_col_t stat_col[];
		public stat_snode_t stat_snode[];
		public int panhows[];
		public cp_panel_t cp_panel[]; /* panels on the critical path */
		public desc_eft_t desc_eft[]; /* all we need to know from descendants */
		public int        cp_firstkid[], cp_nextkid[]; /* linked list of children */
		public int        height[];
		public float      flops_by_height[];
	}

	public static class Branch {
		public int root, first_desc, which_bin;
		public Branch next;
	};


	/* Statistics for supernode and panel size */
	public static int 	no_panels;
	public static float   sum_w;          /* Sum (Wi) */
	public static float 	sum_np_w;       /* Sum (Npi*Wi) */
	public static int 	max_np;
	public static int     no_sups;
	public static float   sum_sup;        /* Sum (Supi) */
	public static int     max_sup;
	public static float reuse_flops;    /* Triangular solve and matrix vector multiply */
	public static float   reuse_data;     /* Doubles in updating supernode */

	/* Statistics for blas operations */
	public static int     num_blas;       /* no of BLAS2 operations, including trsv/gemv */
	public static int     max_blas_n;     /* max dimension n in tri-solve and mat-vec */
	public static int     min_blas_n;     /* min dimension n in tri-solve and mat-vec */
	public static float   sum_blas_n;     /* sum of "        "        " */
	public static int     max_gemv_m;     /* max dimension m in mat-vec */
	public static int     min_gemv_m;     /* max dimension m in mat-vec */
	public static float   sum_gemv_m;     /* sum of "        "        " */
	public static int     lda_blas_m;
	public static int     lda_blas_n;
	public static float gemv_ops[];      /* flops distribution on (m,n) */
	public static float trsv_ops[];      /* flops distribution on n */

	public static double i_trsv_ops(double[] trsv_ops, int i) {
		return trsv_ops[i];
	}
	public static int ij_gemv_ops(int[] gemv_ops, int i, int j) {
		return gemv_ops[j*lda_blas_m + i];
	}

}

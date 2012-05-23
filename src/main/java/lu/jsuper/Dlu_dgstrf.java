/*! @file Dlu_dgstrf.java
 * \brief Computes an LU factorization of a general sparse matrix
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

import lu.jsuper.Dlu_slu_ddefs.GlobalLU_t;
import lu.jsuper.Dlu_slu_util.Dlu_superlu_options_t;
import lu.jsuper.Dlu_slu_util.SuperLUStat_t;
import lu.jsuper.Dlu_superlu_enum_consts.fact_t;
import lu.jsuper.Dlu_supermatrix.NCPformat;
import lu.jsuper.Dlu_supermatrix.SuperMatrix;

import static lu.jsuper.Dlu_sp_ienv.sp_ienv;
import static lu.jsuper.Dlu_slu_util.SUPERLU_MIN;

import static lu.jsuper.Dlu_dmemory.dLUMemInit;
import static lu.jsuper.Dlu_dmemory.dSetRWork;
import static lu.jsuper.Dlu_memory.SetIWork;


public class Dlu_dgstrf {

	private static GlobalLU_t Glu; /* persistent to facilitate multiple factors. */

	/**! \brief
	 *
	 * <pre>
	 * Purpose
	 * =======
	 *
	 * DGSTRF computes an LU factorization of a general sparse m-by-n
	 * matrix A using partial pivoting with row interchanges.
	 * The factorization has the form
	 *     Pr * A = L * U
	 * where Pr is a row permutation matrix, L is lower triangular with unit
	 * diagonal elements (lower trapezoidal if A.nrow > A.ncol), and U is upper
	 * triangular (upper trapezoidal if A.nrow < A.ncol).
	 *
	 * See supermatrix.h for the definition of 'SuperMatrix' structure.
	 *
	 * Arguments
	 * =========
	 *
	 * options (input) superlu_options_t*
	 *         The structure defines the input parameters to control
	 *         how the LU decomposition will be performed.
	 *
	 * A        (input) SuperMatrix*
	 *	    Original matrix A, permuted by columns, of dimension
	 *          (A.nrow, A.ncol). The type of A can be:
	 *          Stype = SLU_NCP; Dtype = SLU_D; Mtype = SLU_GE.
	 *
	 * relax    (input) int
	 *          To control degree of relaxing supernodes. If the number
	 *          of nodes (columns) in a subtree of the elimination tree is less
	 *          than relax, this subtree is considered as one supernode,
	 *          regardless of the row structures of those columns.
	 *
	 * panel_size (input) int
	 *          A panel consists of at most panel_size consecutive columns.
	 *
	 * etree    (input) int*, dimension (A.ncol)
	 *          Elimination tree of A'*A.
	 *          Note: etree is a vector of parent pointers for a forest whose
	 *          vertices are the integers 0 to A.ncol-1; etree[root]==A.ncol.
	 *          On input, the columns of A should be permuted so that the
	 *          etree is in a certain postorder.
	 *
	 * work     (input/output) void*, size (lwork) (in bytes)
	 *          User-supplied work space and space for the output data structures.
	 *          Not referenced if lwork = 0;
	 *
	 * lwork   (input) int
	 *         Specifies the size of work array in bytes.
	 *         = 0:  allocate space internally by system malloc;
	 *         > 0:  use user-supplied work array of length lwork in bytes,
	 *               returns error if space runs out.
	 *         = -1: the routine guesses the amount of space needed without
	 *               performing the factorization, and returns it in
	 *               *info; no other side effects.
	 *
	 * perm_c   (input) int*, dimension (A.ncol)
	 *	    Column permutation vector, which defines the
	 *          permutation matrix Pc; perm_c[i] = j means column i of A is
	 *          in position j in A*Pc.
	 *          When searching for diagonal, perm_c[*] is applied to the
	 *          row subscripts of A, so that diagonal threshold pivoting
	 *          can find the diagonal of A, rather than that of A*Pc.
	 *
	 * perm_r   (input/output) int*, dimension (A.nrow)
	 *          Row permutation vector which defines the permutation matrix Pr,
	 *          perm_r[i] = j means row i of A is in position j in Pr*A.
	 *          If options.Fact = SamePattern_SameRowPerm, the pivoting routine
	 *             will try to use the input perm_r, unless a certain threshold
	 *             criterion is violated. In that case, perm_r is overwritten by
	 *             a new permutation determined by partial pivoting or diagonal
	 *             threshold pivoting.
	 *          Otherwise, perm_r is output argument;
	 *
	 * L        (output) SuperMatrix*
	 *          The factor L from the factorization Pr*A=L*U; use compressed row
	 *          subscripts storage for supernodes, i.e., L has type:
	 *          Stype = SLU_SC, Dtype = SLU_D, Mtype = SLU_TRLU.
	 *
	 * U        (output) SuperMatrix*
	 *	    The factor U from the factorization Pr*A*Pc=L*U. Use column-wise
	 *          storage scheme, i.e., U has types: Stype = SLU_NC,
	 *          Dtype = SLU_D, Mtype = SLU_TRU.
	 *
	 * stat     (output) SuperLUStat_t*
	 *          Record the statistics on runtime and floating-point operation count.
	 *          See slu_util.h for the definition of 'SuperLUStat_t'.
	 *
	 * info     (output) int*
	 *          = 0: successful exit
	 *          < 0: if info = -i, the i-th argument had an illegal value
	 *          > 0: if info = i, and i is
	 *             <= A.ncol: U(i,i) is exactly zero. The factorization has
	 *                been completed, but the factor U is exactly singular,
	 *                and division by zero will occur if it is used to solve a
	 *                system of equations.
	 *             > A.ncol: number of bytes allocated when memory allocation
	 *                failure occurred, plus A.ncol. If lwork = -1, it is
	 *                the estimated amount of space needed, plus A.ncol.
	 *
	 * ======================================================================
	 *
	 * Local Working Arrays:
	 * ======================
	 *   m = number of rows in the matrix
	 *   n = number of columns in the matrix
	 *
	 *   xprune[0:n-1]: xprune[*] points to locations in subscript
	 *	vector lsub[*]. For column i, xprune[i] denotes the point where
	 *	structural pruning begins. I.e. only xlsub[i],..,xprune[i]-1 need
	 *	to be traversed for symbolic factorization.
	 *
	 *   marker[0:3*m-1]: marker[i] = j means that node i has been
	 *	reached when working on column j.
	 *	Storage: relative to original row subscripts
	 *	NOTE: There are 3 of them: marker/marker1 are used for panel dfs,
	 *	      see dpanel_dfs.c; marker2 is used for inner-factorization,
	 *            see dcolumn_dfs.c.
	 *
	 *   parent[0:m-1]: parent vector used during dfs
	 *      Storage: relative to new row subscripts
	 *
	 *   xplore[0:m-1]: xplore[i] gives the location of the next (dfs)
	 *	unexplored neighbor of i in lsub[*]
	 *
	 *   segrep[0:nseg-1]: contains the list of supernodal representatives
	 *	in topological order of the dfs. A supernode representative is the
	 *	last column of a supernode.
	 *      The maximum size of segrep[] is n.
	 *
	 *   repfnz[0:W*m-1]: for a nonzero segment U[*,j] that ends at a
	 *	supernodal representative r, repfnz[r] is the location of the first
	 *	nonzero in this segment.  It is also used during the dfs: repfnz[r]>0
	 *	indicates the supernode r has been explored.
	 *	NOTE: There are W of them, each used for one column of a panel.
	 *
	 *   panel_lsub[0:W*m-1]: temporary for the nonzeros row indices below
	 *      the panel diagonal. These are filled in during dpanel_dfs(), and are
	 *      used later in the inner LU factorization within the panel.
	 *	panel_lsub[]/dense[] pair forms the SPA data structure.
	 *	NOTE: There are W of them.
	 *
	 *   dense[0:W*m-1]: sparse accumulating (SPA) vector for intermediate values;
	 *	    	   NOTE: there are W of them.
	 *
	 *   tempv[0:*]: real temporary used for dense numeric kernels;
	 *	The size of this array is defined by NUM_TEMPV() in slu_ddefs.h.
	 * </pre>
	 */
	public static void dgstrf(Dlu_superlu_options_t options, SuperMatrix A,
	        int relax, int panel_size, int etree[], double work[], int lwork,
	        int perm_c[], int perm_r[], SuperMatrix[] L, SuperMatrix[] U,
	        SuperLUStat_t stat, int[] info) {

	    /* Local working arrays */
	    NCPformat Astore;
	    int       iperm_r[] = null; /* inverse of perm_r; used when
	                                  options.Fact == SamePattern_SameRowPerm */
	    int       iperm_c[]; /* inverse of perm_c */
	    int       iwork[][] = new int[1][];
	    double    dwork[][] = new double[1][];
	    int	      segrep[][] = new int[1][];
	    int       repfnz[][] = new int[1][];
	    int       parent[][] = new int[1][];
	    int       xplore[][] = new int[1][];
	    /* dense[]/panel_lsub[] pair forms a w-wide SPA */
	    int	      panel_lsub[][] = new int[1][];
	    int	      xprune[][] = new int[1][];
	    int	      marker[][] = new int[1][];
	    double    dense[][] = new double[1][];
	    double    tempv[][] = new double[1][];
	    int       relax_end[];
	    double    a[];
	    int       asub[];
	    int       xa_begin[], xa_end[];
	    int       xsup[], supno[];
	    int       xlsub[], xlusup[], xusub[];
	    int       nzlumax;
	    double fill_ratio = sp_ienv(6);  /* estimated fill ratio */

	    /* Local scalars */
	    fact_t    fact = options.Fact;
	    double    diag_pivot_thresh = options.DiagPivotThresh;
	    int       pivrow;   /* pivotal row number in the original matrix A */
	    int       nseg1;	/* no of segments in U-column above panel row jcol */
	    int       nseg;	/* no of segments in each U-column */
	    int jcol;
	    int kcol;	/* end column of a relaxed snode */
	    int icol;
	    int i, k, jj, new_next, iinfo;
	    int       m, n, min_mn, jsupno, fsupc, nextlu, nextu;
	    int       w_def;	/* upper bound on panel width */
	    int       usepr, iperm_r_allocated = 0;
	    int       nnzL, nnzU;
	    int       panel_histo[] = stat.panel_histo;
	    float     ops[] = stat.ops;

	    iinfo    = 0;
	    m        = A.nrow;
	    n        = A.ncol;
	    min_mn   = SUPERLU_MIN(m, n);
	    Astore   = (NCPformat) A.Store;
	    a        = Astore.nzval;
	    asub     = Astore.rowind;
	    xa_begin = Astore.colbeg;
	    xa_end   = Astore.colend;

	    /* Allocate storage common to the factor routines */
	    info[0] = dLUMemInit(fact, work, lwork, m, n, Astore.nnz,
	                       panel_size, fill_ratio, L[0], U[0], Glu, iwork, dwork);
	    if ( info[0] != 0 ) return;

	    xsup    = Glu.xsup;
	    supno   = Glu.supno;
	    xlsub   = Glu.xlsub;
	    xlusup  = Glu.xlusup;
	    xusub   = Glu.xusub;

	    SetIWork(m, n, panel_size, iwork[0], segrep, parent, xplore,
		     repfnz, panel_lsub, xprune, marker);
	    dSetRWork(m, panel_size, dwork[0], dense, tempv);

	}

}

package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.pdgstrf_threadarg_t;
import gov.lbl.superlu.Dlu_pdsp_defs.pxgstrf_shared_t;
import gov.lbl.superlu.Dlu_slu_mt_util.Gstat_t;
import gov.lbl.superlu.Dlu_slu_mt_util.superlumt_options_t;
import gov.lbl.superlu.Dlu_supermatrix.SuperMatrix;

import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_ABORT;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.FACT;

import static gov.lbl.superlu.Dlu.printf;
import static gov.lbl.superlu.Dlu.fprintf;
import static gov.lbl.superlu.Dlu.stderr;
import static gov.lbl.superlu.Dlu.PRNTlevel;

import static gov.lbl.superlu.Dlu_superlu_timer.SuperLU_timer_;
import static gov.lbl.superlu.Dlu_superlu_timer.usertimer_;

import static gov.lbl.superlu.Dlu_pdgstrf_thread_init.pdgstrf_thread_init;
import static gov.lbl.superlu.Dlu_pdgstrf_thread_finalize.pdgstrf_thread_finalize;

import static gov.lbl.superlu.Dlu_pdgstrf_thread.pdgstrf_thread;


public class Dlu_pdgstrf {

	static
	void
	pdgstrf(superlumt_options_t superlumt_options, SuperMatrix A, int perm_r[],
		SuperMatrix L, SuperMatrix U, Gstat_t Gstat, int info[])
	{
	/*
	 * -- SuperLU MT routine (version 2.0) --
	 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
	 * and Xerox Palo Alto Research Center.
	 * September 10, 2007
	 *
	 * Purpose
	 * =======
	 *
	 * PDGSTRF computes an LU factorization of a general sparse nrow-by-ncol
	 * matrix A using partial pivoting with row interchanges. The factorization
	 * has the form
	 *     Pr * A = L * U
	 * where Pr is a row permutation matrix, L is lower triangular with unit
	 * diagonal elements (lower trapezoidal if A.nrow > A.ncol), and U is
	 * upper triangular (upper trapezoidal if A.nrow < A.ncol).
	 *
	 * Arguments
	 * =========
	 *
	 * superlumt_options (input) superlumt_options_t*
	 *        The structure defines the parameters to control how the sparse
	 *        LU factorization is performed. The following fields must be set
	 *        by the user:
	 *
	 *        o nprocs (int)
	 *          Number of processes to be spawned and used for factorization.
	 *
	 *        o refact (yes_no_t)
	 *          Specifies whether this is first time or subsequent factorization.
	 *          = NO:  this factorization is treated as the first one;
	 *          = YES: it means that a factorization was performed prior to this
	 *                 one. Therefore, this factorization will re-use some
	 *                 existing data structures, such as L and U storage, column
	 *                 elimination tree, and the symbolic information of the
	 *                 Householder matrix.
	 *
	 *        o panel_size (int)
	 *          A panel consists of at most panel_size consecutive columns.
	 *
	 *        o relax (int)
	 *          Degree of relaxing supernodes. If the number of nodes (columns)
	 *          in a subtree of the elimination tree is less than relax, this
	 *          subtree is considered as one supernode, regardless of the row
	 *          structures of those columns.
	 *
	 *        o diag_pivot_thresh (double)
	 *	    Diagonal pivoting threshold. At step j of Gaussian elimination,
	 *          if abs(A_jj) >= diag_pivot_thresh * (max_(i>=j) abs(A_ij)),
	 *          use A_jj as pivot. 0 <= diag_pivot_thresh <= 1. The default
	 *          value is 1.0, corresponding to partial pivoting.
	 *
	 *        o usepr (yes_no_t)
	 *          Whether the pivoting will use perm_r specified by the user.
	 *          = YES: use perm_r; perm_r is input, unchanged on exit.
	 *          = NO:  perm_r is determined by partial pivoting, and is output.
	 *
	 *        o drop_tol (double) (NOT IMPLEMENTED)
	 *	    Drop tolerance parameter. At step j of the Gaussian elimination,
	 *          if abs(A_ij)/(max_i abs(A_ij)) < drop_tol, drop entry A_ij.
	 *          0 <= drop_tol <= 1. The default value of drop_tol is 0,
	 *          corresponding to not dropping any entry.
	 *
	 *        o perm_c (int*)
	 *	    Column permutation vector of size A.ncol, which defines the
	 *          permutation matrix Pc; perm_c[i] = j means column i of A is
	 *          in position j in A*Pc.
	 *
	 *        o perm_r (int*)
	 *	    Column permutation vector of size A.nrow.
	 *          If superlumt_options.usepr = NO, this is an output argument.
	 *
	 *        o work (void*) of size lwork
	 *          User-supplied work space and space for the output data structures.
	 *          Not referenced if lwork = 0;
	 *
	 *        o lwork (int)
	 *          Specifies the length of work array.
	 *            = 0:  allocate space internally by system malloc;
	 *            > 0:  use user-supplied work array of length lwork in bytes,
	 *                  returns error if space runs out.
	 *            = -1: the routine guesses the amount of space needed without
	 *                  performing the factorization, and returns it in
	 *                  superlu_memusage.total_needed; no other side effects.
	 *
	 * A      (input) SuperMatrix*
	 *	  Original matrix A, permuted by columns, of dimension
	 *        (A.nrow, A.ncol). The type of A can be:
	 *        Stype = NCP; Dtype = _D; Mtype = GE.
	 *
	 * perm_r (input/output) int*, dimension A.nrow
	 *        Row permutation vector which defines the permutation matrix Pr,
	 *        perm_r[i] = j means row i of A is in position j in Pr*A.
	 *        If superlumt_options.usepr = NO, perm_r is output argument;
	 *        If superlumt_options.usepr = YES, the pivoting routine will try
	 *           to use the input perm_r, unless a certain threshold criterion
	 *           is violated. In that case, perm_r is overwritten by a new
	 *           permutation determined by partial pivoting or diagonal
	 *           threshold pivoting.
	 *
	 * L      (output) SuperMatrix*
	 *        The factor L from the factorization Pr*A=L*U; use compressed row
	 *        subscripts storage for supernodes, i.e., L has type:
	 *        Stype = SCP, Dtype = _D, Mtype = TRLU.
	 *
	 * U      (output) SuperMatrix*
	 *	  The factor U from the factorization Pr*A*Pc=L*U. Use column-wise
	 *        storage scheme, i.e., U has types: Stype = NCP, Dtype = _D,
	 *        Mtype = TRU.
	 *
	 * Gstat  (output) Gstat_t*
	 *        Record all the statistics about the factorization;
	 *        See Gstat_t structure defined in slu_mt_util.h.
	 *
	 * info   (output) int*
	 *        = 0: successful exit
	 *        < 0: if info = -i, the i-th argument had an illegal value
	 *        > 0: if info = i, and i is
	 *             <= A.ncol: U(i,i) is exactly zero. The factorization has
	 *                been completed, but the factor U is exactly singular,
	 *                and division by zero will occur if it is used to solve a
	 *                system of equations.
	 *             > A.ncol: number of bytes allocated when memory allocation
	 *                failure occurred, plus A.ncol.
	 *
	 */
	    pdgstrf_threadarg_t pdgstrf_threadarg[];
	    pxgstrf_shared_t pxgstrf_shared = new pxgstrf_shared_t();
	    int nprocs = superlumt_options.nprocs;
	    int i;
		int iinfo;
	    double    utime[] = Gstat.utime;
	    double    usrtime, wtime;
	    Thread thread_id[];
	    Object      status[];
//	    void      *pdgstrf_thread(void *);


	    /* --------------------------------------------------------------
	       Initializes the parallel data structures for pdgstrf_thread().
	       --------------------------------------------------------------*/
	    pdgstrf_threadarg = pdgstrf_thread_init(A, L, U, superlumt_options,
						    pxgstrf_shared, Gstat, info);
	    if ( info[0] != 0 ) return;

	    /* Start timing factorization. */
	    usrtime = usertimer_();
	    wtime = SuperLU_timer_();

	    /* ------------------------------------------------------------
	       Use POSIX threads.
	       ------------------------------------------------------------*/
	if ( true )	{  /* Use pthread ... */

	    /* Create nproc threads for concurrent factorization. */
	    thread_id = new Thread [nprocs];

	    for (final pdgstrf_threadarg_t arg : pdgstrf_threadarg) {
	    	new Thread() {
				public void run() {
					pdgstrf_thread(arg);
				}
			}.start();
		}

//	    for (i = 0; i < nprocs; ++i) {
//		if ( iinfo = pthread_create(&thread_id[i],
//					    null,
//					    pdgstrf_thread,
//					    &(pdgstrf_threadarg[i])) ) {
//		    fprintf(stderr, "pthread_create: %d\n", iinfo);
//		    SUPERLU_ABORT("pthread_create()");
//		}
//	    }
//
//	    /* Wait for all threads to terminate. */
//	    for (i = 0; i < nprocs; i++)
//		pthread_join(thread_id[i], &status);
	/* _PTHREAD */

	    /* ------------------------------------------------------------
	       On all other systems, use single processor.
	       ------------------------------------------------------------*/
	} else {

	    printf("pdgstrf() is not parallelized on this machine.\n");
	    printf("pdgstrf() will be run on single processor.\n");
	    pdgstrf_thread( pdgstrf_threadarg[0] );

	}

	    wtime = SuperLU_timer_() - wtime;
	    usrtime = usertimer_() - usrtime;
	    utime[FACT.ordinal()] = wtime;

	if ( PRNTlevel==1 ) {
	    printf(".. pdgstrf_thread() returns info %d, usrtime %.2f, wtime %.2f\n",
	           info[0], usrtime, wtime);
	}

	    /* check_mem_leak("after pdgstrf_thread()"); */

	    /* ------------------------------------------------------------
	       Clean up and free storage after multithreaded factorization.
	       ------------------------------------------------------------*/
	    pdgstrf_thread_finalize(pdgstrf_threadarg, pxgstrf_shared,
				    A, perm_r, L, U);

	}

}

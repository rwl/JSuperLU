package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_pdsp_defs.pdgstrf_threadarg_t;
import gov.lbl.superlu.Dlu_pdsp_defs.pxgstrf_shared_t;
import gov.lbl.superlu.Dlu_slu_mt_util.superlumt_options_t;
import gov.lbl.superlu.Dlu_supermatrix.NCPformat;
import gov.lbl.superlu.Dlu_supermatrix.SCPformat;
import gov.lbl.superlu.Dlu_supermatrix.SuperMatrix;

import static gov.lbl.superlu.Dlu_util.countnz;
import static gov.lbl.superlu.Dlu_util.fixupL;
import static gov.lbl.superlu.Dlu_util.compressSUP;
import static gov.lbl.superlu.Dlu_util.PrintGLGU;
import static gov.lbl.superlu.Dlu_util.PrintInt10;

import static gov.lbl.superlu.Dlu.COMPRESS_LUSUP;

import static gov.lbl.superlu.Dlu_slu_mt_util.yes_no_t.YES;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_MIN;

import static gov.lbl.superlu.Dlu_supermatrix.Stype_t.SLU_SCP;
import static gov.lbl.superlu.Dlu_supermatrix.Dtype_t.SLU_D;
import static gov.lbl.superlu.Dlu_supermatrix.Mtype_t.SLU_TRLU;
import static gov.lbl.superlu.Dlu_supermatrix.Stype_t.SLU_NCP;
import static gov.lbl.superlu.Dlu_supermatrix.Mtype_t.SLU_TRU;

import static gov.lbl.superlu.Dlu_pdutil.dCreate_SuperNode_Permuted;
import static gov.lbl.superlu.Dlu_pdutil.dCreate_CompCol_Permuted;

import static gov.lbl.superlu.Dlu.printf;
import static gov.lbl.superlu.Dlu.DEBUGlevel;

import static gov.lbl.superlu.Dlu_pxgstrf_synch.ParallelFinalize;
import static gov.lbl.superlu.Dlu_pxgstrf_synch.QueryQueue;


public class Dlu_pdgstrf_thread_finalize {

	static
	void
	pdgstrf_thread_finalize(pdgstrf_threadarg_t pdgstrf_threadarg[],
				pxgstrf_shared_t pxgstrf_shared,
				SuperMatrix A, int perm_r[],
				SuperMatrix L, SuperMatrix U
				)
	{
	/*
	 * -- SuperLU MT routine (version 2.0) --
	 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
	 * and Xerox Palo Alto Research Center.
	 * September 10, 2007
	 *
	 *
	 * Purpose
	 * =======
	 *
	 * pdgstrf_thread_finalize() performs cleanups after the multithreaded
	 * factorization pdgstrf_thread(). It sets up the L and U data
	 * structures, and deallocats the storage associated with the structures
	 * pxgstrf_shared and pdgstrf_threadarg.
	 *
	 * Arguments
	 * =========
	 *
	 * pdgstrf_threadarg (input) pdgstrf_threadarg_t*
	 *          The structure contains the parameters to each thread.
	 *
	 * pxgstrf_shared (input) pxgstrf_shared_t*
	 *          The structure contains the shared task queue, the
	 *          synchronization variables, and the L and U data structures.
	 *
	 * A        (input) SuperMatrix*
	 *	    Original matrix A, permutated by columns, of dimension
	 *          (A.nrow, A.ncol). The type of A can be:
	 *          Stype = NCP; Dtype = _D; Mtype = GE.
	 *
	 * perm_r   (input) int*, dimension A.nrow
	 *          Row permutation vector which defines the permutation matrix Pr,
	 *          perm_r[i] = j means row i of A is in position j in Pr*A.
	 *
	 * L        (output) SuperMatrix*
	 *          The factor L from the factorization Pr*A=L*U; use compressed row
	 *          subscripts storage for supernodes, i.e., L has type:
	 *          Stype = SCP, Dtype = _D, Mtype = TRLU.
	 *
	 * U        (output) SuperMatrix*
	 *	    The factor U from the factorization Pr*A*Pc=L*U. Use column-wise
	 *          storage scheme, i.e., U has type:
	 *          Stype = NCP, Dtype = _D, Mtype = TRU.
	 *
	 *
	 */
	    int nprocs, n, i, iinfo;
	    int[]       nnzL, nnzU;
	    nnzL = new int[1];
	    nnzU = new int[1];
	    superlumt_options_t superlumt_options;
	    GlobalLU_t Glu;

	    n = A.ncol;
	    superlumt_options = pdgstrf_threadarg[0].superlumt_options;
	    Glu = pxgstrf_shared.Glu;
	    Glu.supno[n] = Glu.nsuper;

	    countnz(n, pxgstrf_shared.xprune, nnzL, nnzU, Glu);
	    fixupL(n, perm_r, Glu);

	if (COMPRESS_LUSUP) {
	    compressSUP(n, pxgstrf_shared.Glu);
	}

	    if ( superlumt_options.refact == YES ) {
	        /* L and U structures may have changed due to possibly different
		   pivoting, although the storage is available. */
	        ((SCPformat)L.Store).nnz = nnzL[0];
		((SCPformat)L.Store).nsuper = Glu.supno[n];
		((NCPformat)U.Store).nnz = nnzU[0];
	    } else {
		dCreate_SuperNode_Permuted(L, A.nrow, A.ncol, nnzL[0], Glu.lusup,
					   Glu.xlusup, Glu.xlusup_end,
					   Glu.lsub, Glu.xlsub, Glu.xlsub_end,
					   Glu.supno, Glu.xsup, Glu.xsup_end,
					   SLU_SCP, SLU_D, SLU_TRLU);
		dCreate_CompCol_Permuted(U, A.nrow, A.ncol, nnzU[0], Glu.ucol,
					 Glu.usub, Glu.xusub, Glu.xusub_end,
					 SLU_NCP, SLU_D, SLU_TRU);
	    }

	    /* Combine the INFO returned from individual threads. */
	    iinfo = 0;
	    nprocs = superlumt_options.nprocs;
	    for (i = 0; i < nprocs; ++i) {
	        if ( pdgstrf_threadarg[i].info[0] != 0 ) {
		    if (iinfo != 0) iinfo=SUPERLU_MIN(iinfo, pdgstrf_threadarg[i].info[0]);
		    else iinfo = pdgstrf_threadarg[i].info[0];
		}
	    }
	    pxgstrf_shared.info = iinfo;

	if ( DEBUGlevel>=2 ) {
	    printf("Last nsuper %d\n", Glu.nsuper);
	    QueryQueue(pxgstrf_shared.taskq);
	    PrintGLGU(n, pxgstrf_shared.xprune, Glu);
	    PrintInt10("perm_r", n, perm_r);
	    PrintInt10("inv_perm_r", n, pxgstrf_shared.inv_perm_r);
	}

	    /* Deallocate the storage used by the parallel scheduling algorithm. */
	    ParallelFinalize(pxgstrf_shared);
	    pdgstrf_threadarg = null;
	    pxgstrf_shared.inv_perm_r = null;
	    pxgstrf_shared.inv_perm_c = null;
	    pxgstrf_shared.xprune = null;
	    pxgstrf_shared.ispruned = null;

	if ( DEBUGlevel>=1 ) {
	    printf("** pdgstrf_thread_finalize() called\n");
	}
	}

}

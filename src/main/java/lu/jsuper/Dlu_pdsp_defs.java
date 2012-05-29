/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 * Sparse matrix types and function prototypes.
 *
 */
package lu.jsuper;

import lu.jsuper.Dlu_pxgstrf_synch.pan_status_t;
import lu.jsuper.Dlu_slu_mt_util.Gstat_t;
import lu.jsuper.Dlu_slu_mt_util.superlumt_options_t;
import lu.jsuper.Dlu_supermatrix.SuperMatrix;

public class Dlu_pdsp_defs {

	/*
	 * *************************************************
	 *  Global data structures used in LU factorization
	 * *************************************************
	 *
	 *   nsuper: number of supernodes = nsuper+1, numbered between 0 and nsuper.
	 *
	 *   (supno, xsup, xsup_end):
	 *      supno[i] is the supernode number to which column i belongs;
	 *	xsup[s] points to the first column of supernode s;
	 *      xsup_end[s] points to one past the last column of supernode s.
	 *	Example: supno  0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
	 *	          xsup  0 1 2 4 7
	 *            xsup_end  1 2 4 7 12
	 *	Note: dfs will be performed on supernode rep. relative to the new
	 *	      row pivoting ordering
	 *
	 *   (lsub, xlsub, xlsub_end):
	 *      lsub[*] contains the compressed subscripts of the supernodes;
	 *      xlsub[j] points to the starting location of the j-th column in
	 *               lsub[*];
	 *      xlsub_end[j] points to one past the ending location of the j-th
	 *               column in lsub[*].
	 *	Storage: original row subscripts in A.
	 *
	 *      During the course of sparse LU factorization, we also use
	 *	(lsub, xlsub, xlsub_end, xprune) to represent symmetrically
	 *      pruned graph. Contention will occur when one processor is
	 *      performing DFS on supernode S, while another processor is pruning
	 *      supernode S. We use the following data structure to deal with
	 *      this problem. Suppose each supernode contains columns {s,s+1,...,t},
	 *      with first column s and last column t.
	 *
	 *      (1) if t > s, only the subscript sets for column s and column t
	 *          are stored. Column t represents pruned adjacency structure.
	 *
	 *                  --------------------------------------------
	 *          lsub[*]    ... |   col s    |   col t   | ...
	 *                  --------------------------------------------
	 *                          ^            ^           ^
	 *                       xlsub[s]    xlsub_end[s]  xlsub_end[s+1]
	 *                                   xlsub[s+1]      :
	 *                                       :           :
	 *                                       :         xlsub_end[t]
	 *                                   xlsub[t]      xprune[t]
	 *                                   xprune[s]
	 *
	 *      (2) if t == s, i.e., a singleton supernode, the subscript set
	 *          is stored twice:
	 *
	 *                  --------------------------------------
	 *          lsub[*]    ... |      s     |     s     | ...
	 *                  --------------------------------------
	 *                          ^            ^           ^
	 *                       xlsub[s]   xlsub_end[s]  xprune[s]
	 *
	 *      There are two subscript sets for each supernode, the last column
	 *      structures (for pruning) will be removed after the numerical LU
	 *      factorization phase:
	 *        o lsub[j], j = xlsub[s], ..., xlsub_end[s]-1
	 *          is the structure of column s (i.e. structure of this supernode).
	 *          It is used for the storage of numerical values.
	 *	  o lsub[j], j = xlsub[t], ..., xlsub_end[t]-1
	 *	    is the structure of the last column t of this supernode.
	 *	    It is for the purpose of symmetric pruning. Therefore, the
	 *	    structural subscripts can be rearranged without making physical
	 *	    interchanges among the numerical values.
	 *
	 *       DFS will traverse the first subscript set if the supernode
	 *       has not been pruned; otherwise it will traverse the second
	 *       subscript list, i.e., the part of the pruned graph.
	 *
	 *   (lusup, xlusup, xlusup_end):
	 *      lusup[*] contains the numerical values of the supernodes;
	 *      xlusup[j] points to the starting location of the j-th column in
	 *                storage vector lusup[*];
	 *      xlusup_end[j] points to one past the ending location of the j-th
	 *                column in lusup[*].
	 *	Each supernode is stored in column-major, consistent with Fortran
	 *      two-dimensional array storage.
	 *
	 *   (ucol, usub, xusub, xusub_end):
	 *      ucol[*] stores the numerical values of the U-columns above the
	 *              supernodes.
	 *      usub[k] stores the row subscripts of nonzeros ucol[k];
	 *      xusub[j] points to the starting location of column j in ucol/usub[];
	 *      xusub_end[j] points to one past the ending location column j in
	 *                   ucol/usub[].
	 *	Storage: new row subscripts; that is indexed intp PA.
	 *
	 */
	static class GlobalLU_t {
	    int     xsup[];    /* supernode and column mapping */
	    int     xsup_end[];
	    int     supno[];
	    int     lsub[];    /* compressed L subscripts */
	    int	    xlsub[];
	    int     xlsub_end[];
	    double  lusup[];   /* L supernodes */
	    int     xlusup[];
	    int     xlusup_end[];
	    double  ucol[];    /* U columns */
	    int     usub[];
	    int	    xusub[];
	    int     xusub_end[];
	    int     nsuper;   /* current supernode number */
	    int     nextl;    /* next position in lsub[] */
	    int     nextu;    /* next position in usub[]/ucol[] */
	    int     nextlu;   /* next position in lusup[] */
	    int     nzlmax;   /* current max size of lsub[] */
	    int     nzumax;   /*    "    "    "      ucol[] */
	    int     nzlumax;  /*    "    "    "     lusup[] */
	    /* ---------------------------------------------------------------
	     *  Memory managemant for L supernodes
	     */
	    int  map_in_sup[];  /* size n+1 - the address offset of each column
	                         * in lusup[*], which is divided into regions
				* by the supernodes of Householder matrix H.
				* If column k starts a supernode in H,
				* map_in_sup[k] is the next open position in
				* lusup[*]; otherwise map_in_sup[k] gives the
				* offset (negative) to the leading column
				* of the supernode in H.
				*/
	    int  dynamic_snode_bound;
	    /* --------------------------------------------------------------- */
	}


	/*
	 * *********************************************************************
	 * The pxgstrf_shared_t structure contains the shared task queue and
	 * the synchronization variables to facilitate parallel factorization.
	 * It also contains the shared L and U data structures.
	 * *********************************************************************
	 */
	static class pxgstrf_shared_t {
	    /* ----------------------------------------------------------------
	     * Global variables introduced in parallel code for synchronization.
	     */
	    volatile int tasks_remain; /* number of untaken panels */
	    int          num_splits;   /* number of panels split at the top */
	    queue_t      taskq;        /* size ncol - shared work queue */
	    mutex_t      lu_locks[];    /* 5 named mutual exclusive locks */
	    volatile int spin_locks[];  /* size ncol - mark every busy column */
	    pan_status_t pan_status[];  /* size ncol - panel status */
	    int          fb_cols[];     /* size ncol - mark farthest busy column */
	    /* ---------------------------------------------------------------- */
	    int        inv_perm_c[];
	    int        inv_perm_r[];
	    int        xprune[];
	    int        ispruned[];
	    SuperMatrix A;
	    GlobalLU_t Glu;
	    Gstat_t    Gstat;
	    int        info;
	}

	/* Arguments passed to each thread. */
	static class pdgstrf_threadarg_t {
	    int  pnum; /* process number */
	    int  info; /* error code returned from each thread */
	    superlumt_options_t superlumt_options;
	    pxgstrf_shared_t  pxgstrf_shared; /* shared for LU factorization */
	}

}

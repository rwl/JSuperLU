/*! @file Dlu_slu_ddefs.java
 * \brief Header file for real operations
 *
 * <pre>
 * -- SuperLU routine (version 4.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November, 2010
 *
 * Global data structures used in LU factorization -
 *
 *   nsuper: #supernodes = nsuper + 1, numbered [0, nsuper].
 *   (xsup,supno): supno[i] is the supernode no to which i belongs;
 *	xsup(s) points to the beginning of the s-th supernode.
 *	e.g.   supno 0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	        xsup 0 1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new
 *	      row pivoting ordering
 *
 *   (xlsub,lsub): lsub[*] contains the compressed subscript of
 *	rectangular supernodes; xlsub[j] points to the starting
 *	location of the j-th column in lsub[*]. Note that xlsub
 *	is indexed by column.
 *	Storage: original row subscripts
 *
 *      During the course of sparse LU factorization, we also use
 *	(xlsub,lsub) for the purpose of symmetric pruning. For each
 *	supernode {s,s+1,...,t=s+r} with first column s and last
 *	column t, the subscript set
 *		lsub[j], j=xlsub[s], .., xlsub[s+1]-1
 *	is the structure of column s (i.e. structure of this supernode).
 *	It is used for the storage of numerical values.
 *	Furthermore,
 *		lsub[j], j=xlsub[t], .., xlsub[t+1]-1
 *	is the structure of the last column t of this supernode.
 *	It is for the purpose of symmetric pruning. Therefore, the
 *	structural subscripts can be rearranged without making physical
 *	interchanges among the numerical values.
 *
 *	However, if the supernode has only one column, then we
 *	only keep one set of subscripts. For any subscript interchange
 *	performed, similar interchange must be done on the numerical
 *	values.
 *
 *	The last column structures (for pruning) will be removed
 *	after the numercial LU factorization phase.
 *
 *   (xlusup,lusup): lusup[*] contains the numerical values of the
 *	rectangular supernodes; xlusup[j] points to the starting
 *	location of the j-th column in storage vector lusup[*]
 *	Note: xlusup is indexed by column.
 *	Each rectangular supernode is stored by column-major
 *	scheme, consistent with Fortran 2-dim array storage.
 *
 *   (xusub,ucol,usub): ucol[*] stores the numerical values of
 *	U-columns outside the rectangular supernodes. The row
 *	subscript of nonzero ucol[k] is stored in usub[k].
 *	xusub[i] points to the starting location of column i in ucol.
 *	Storage: new row subscripts; that is subscripts of PA.
 * </pre>
 */
package lu.jsuper;

import lu.jsuper.Dlu_slu_util.ExpHeader;
import lu.jsuper.Dlu_slu_util.LU_stack_t;
import lu.jsuper.Dlu_superlu_enum_consts.LU_space_t;


public class Dlu_slu_ddefs {

	public static class GlobalLU_t {
	    int     xsup[];    /* supernode and column mapping */
	    int     supno[];
	    int     lsub[];    /* compressed L subscripts */
	    int	    xlsub[];
	    double  lusup[];   /* L supernodes */
	    int     xlusup[];
	    double  ucol[];    /* U columns */
	    int     usub[];
	    int	    xusub[];
	    int     nzlmax;   /* current max size of lsub */
	    int     nzumax;   /*    "    "    "      ucol */
	    int     nzlumax;  /*    "    "    "     lusup */
	    int     n;        /* number of columns in the matrix */
	    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
	    int     num_expansions;
	    ExpHeader expanders[]; /* Array of pointers to 4 types of memory */
	    LU_stack_t stack;     /* use user supplied memory */
	}

}

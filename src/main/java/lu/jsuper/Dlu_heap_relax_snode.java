/*
 * -- SuperLU_MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 */
package lu.jsuper;

import lu.jsuper.Dlu_pxgstrf_synch.pxgstrf_relax_t;
import lu.jsuper.Dlu_slu_mt_util.superlumt_options_t;

import static lu.jsuper.Dlu_slu_mt_util.SUPERLU_ABORT;
import static lu.jsuper.Dlu_slu_mt_util.SUPERLU_MIN;
import static lu.jsuper.Dlu_sp_coletree.TreePostorder;

import static lu.jsuper.Dlu.DOMAINS;
import static lu.jsuper.Dlu.printf;
import static lu.jsuper.Dlu.PRNTlevel;


public class Dlu_heap_relax_snode {

	static
	void
	heap_relax_snode (
			  final     int n,
			  superlumt_options_t superlumt_options,
			  pxgstrf_relax_t pxgstrf_relax[] /* relaxed s-nodes */
			  )
	{
	/*
	 * Purpose
	 * =======
	 *    heap_relax_snode() - Identify the initial relaxed supernodes, assuming
	 *    that the matrix has been reordered according to the postorder of
	 *    the etree.
	 *
	 */
	    int i, j, k, l, parent;
	    int snode_start;	/* beginning of a snode */
	    int et_save[], post[], inv_post[], iwork[];
	    int nsuper_et = 0, nsuper_et_post = 0;

	    int fcol;	 /* beginning of a snode */
	    int desc[];  /* no of descendants of each etree node. */
	    int et[] = superlumt_options.etree; /* column elimination tree */
	    int relax = superlumt_options.relax; /* maximum no of columns allowed
						     in a relaxed s-node */

	    desc = intCalloc(n+1);

	    /* The etree may not be postordered, but is always heap-ordered. */

	    if ( (iwork = (int[]) intMalloc(3*n+2)) == null )
		SUPERLU_ABORT("SUPERLU_MALLOC fails for iwork[]");
	    inv_post = iwork    + n+1;
	    et_save  = inv_post + n+1;

	    /* Post order etree */
	    post = (int []) TreePostorder(n, et);
	    for (i = 0; i < n+1; ++i) inv_post[post[i]] = i;

	    /* Renumber etree in postorder */
	    for (i = 0; i < n; ++i) {
	        iwork[post[i]] = post[et[i]];
		et_save[i] = et[i]; /* Save the original etree */
	    }
	    for (i = 0; i < n; ++i) et[i] = iwork[i];

	    /* Compute the number of descendants of each node in the etree */
	    for (j = 0; j < n; j++) desc[j] = 0;
	    for (j = 0; j < n; j++) {
		parent = et[j];
		if ( parent != n )  /* not the dummy root */
		    desc[parent] += desc[j] + 1;
	    }

	    /* Identify the relaxed supernodes by postorder traversal of the etree. */
	    for (j = 0; j < n; ) {
	     	parent = et[j];
	        snode_start = j;
	 	while ( parent != n && desc[parent] < relax ) {
		    j = parent;
		    parent = et[j];
		}
		/* Found a supernode in postordered etree; j is the last column. */
		++nsuper_et_post;
		k = n;
		for (i = snode_start; i <= j; ++i) k = SUPERLU_MIN(k, inv_post[i]);
		l = inv_post[j];
		if ( (l - k) == (j - snode_start) ) {
		    /* It's also a supernode in the original etree */
		    pxgstrf_relax[nsuper_et].fcol = snode_start;
		    pxgstrf_relax[nsuper_et].size = j - snode_start + 1;
	if (DOMAINS) {
			throw new UnsupportedOperationException();
//		    for (i = snode_start; i <= j; ++i) in_domain[i] = RELAXED_SNODE;
	}
		    ++nsuper_et;
		} else {
		    for (i = snode_start; i <= j; ++i) {
		        l = inv_post[i];
		        if ( desc[i] == 0 ) { /* leaf */
			    pxgstrf_relax[nsuper_et].fcol = l; /* relax_end[l] = l;*/
			    pxgstrf_relax[nsuper_et].size = 1;
			    ++nsuper_et;
			}
		    }
		}
		j++;
		/* Search for a new leaf */
		while ( desc[j] != 0 && j < n ) j++;
	    }

	if ( PRNTlevel>=1 ) {
	    printf(".. heap_snode_relax:\n" +
		   "\tNo of relaxed snodes in postordered etree:\t%d\n" +
		   "\tNo of relaxed snodes in original etree:\t%d\n",
		   nsuper_et_post, nsuper_et);
	}

	    /* Recover the original etree */
	    for (i = 0; i < n; ++i) et[i] = et_save[i];
	}

}

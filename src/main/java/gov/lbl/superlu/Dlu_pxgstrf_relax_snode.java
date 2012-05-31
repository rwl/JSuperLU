package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pxgstrf_synch.pxgstrf_relax_t;
import gov.lbl.superlu.Dlu_slu_mt_util.superlumt_options_t;

import static gov.lbl.superlu.Dlu.DOMAINS;
import static gov.lbl.superlu.Dlu.PRNTlevel;
import static gov.lbl.superlu.Dlu.printf;

import static gov.lbl.superlu.Dlu_pmemory.intCalloc;


public class Dlu_pxgstrf_relax_snode {

	static
	void
	pxgstrf_relax_snode(
			    final int n, /* number of columns in the matrix */
			    superlumt_options_t superlumt_options,
			    pxgstrf_relax_t pxgstrf_relax[] /* relaxed s-nodes */
			    )
	{
	/*
	 * -- SuperLU MT routine (version 2.0) --
	 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
	 * and Xerox Palo Alto Research Center.
	 * September 10, 2007
	 *
	 * Purpose
	 * =======
	 *   pxgstrf_relax_snode() identifes the initial relaxed supernodes,
	 *   assuming that the matrix has been reordered according to the postorder
	 *   of the etree.
	 *
	 */
	    int j, parent, rs;
	    int fcol;	 /* beginning of a snode */
	    int desc[];  /* no of descendants of each etree node. */
	    int etree[] = superlumt_options.etree; /* column elimination tree */
	    int relax = superlumt_options.relax; /* maximum no of columns allowed
						     in a relaxed s-node */

	    desc = intCalloc(n+1);

	    /* Compute the number of descendants of each node in the etree */
	    for (j = 0; j < n; j++) {
		parent = etree[j];
		desc[parent] += desc[j] + 1;
	    }

	    rs = 1;

	    /* Identify the relaxed supernodes by postorder traversal of the etree. */
	    for (j = 0; j < n; ) {
	     	parent = etree[j];
	        fcol = j;
	 	while ( parent != n && desc[parent] < relax ) {
		    j = parent;
		    parent = etree[j];
		}
		/* found a supernode with j being the last column. */
		pxgstrf_relax[rs].fcol = fcol;
		pxgstrf_relax[rs].size = j - fcol + 1;
	if (DOMAINS) {
		throw new UnsupportedOperationException();
//		for (i = fcol; i <= j; ++i) in_domain[i] = RELAXED_SNODE;
	}
		j++;    rs++;
		/* Search for a new leaf */
		while ( desc[j] != 0 && j < n ) j++;
	    }

	    pxgstrf_relax[rs].fcol = n;
	    pxgstrf_relax[0].size = rs-1; /* number of relaxed supernodes */

	if (PRNTlevel==1) {
	    printf(".. No of relaxed s-nodes %d\n", pxgstrf_relax[0].size);
	}

	}

}

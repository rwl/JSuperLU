package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_slu_mt_util.superlumt_options_t;
import gov.lbl.superlu.Dlu_supermatrix.SuperMatrix;

import static gov.lbl.superlu.Dlu.DEBUGlevel;
import static gov.lbl.superlu.Dlu.printf;
import static gov.lbl.superlu.Dlu_util.Destroy_CompCol_Permuted;



public class Dlu_pxgstrf_finalize {

	static
	void
	pxgstrf_finalize(superlumt_options_t superlumt_options, SuperMatrix AC)
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
	 * pxgstrf_finalize() deallocates storage after factorization pxgstrf().
	 *
	 * Arguments
	 * =========
	 *
	 * superlumt_options (input) superlumt_options_t*
	 *        The structure contains the parameters to facilitate sparse
	 *        LU factorization.
	 *
	 * AC (input) SuperMatrix*
	 *        The original matrix with columns permuted.
	 */
	    superlumt_options.etree = null;
	    superlumt_options.colcnt_h = null;
	    superlumt_options.part_super_h = null;
	    Destroy_CompCol_Permuted(AC);
	if ( DEBUGlevel>=1 ) {
	    printf("** pxgstrf_finalize() called\n");
	}
	}

}

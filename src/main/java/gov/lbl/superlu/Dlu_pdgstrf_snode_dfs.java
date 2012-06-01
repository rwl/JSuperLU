package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_pdsp_defs.pxgstrf_shared_t;

import static gov.lbl.superlu.Dlu_pmemory.Glu_alloc;
import static gov.lbl.superlu.Dlu_pxgstrf_synch.NewNsuper;
import static gov.lbl.superlu.Dlu_slu_mt_util.MemType.LSUB;


public class Dlu_pdgstrf_snode_dfs {

	static
	int
	pdgstrf_snode_dfs(
			  final int  pnum,      /* process number */
			  final int  jcol,	  /* in - start of the supernode */
			  final int  kcol, 	  /* in - end of the supernode */
			  final int  asub[],     /* in */
			  final int  xa_begin[], /* in */
			  final int  xa_end[],   /* in */
			  int        xprune[],   /* out */
			  int        marker[],   /* modified */
			  int        col_lsub[], /* values are irrelevant on entry
						   and on return */
			  pxgstrf_shared_t pxgstrf_shared /* modified */
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
	 *    pdgstrf_snode_dfs() determines the union of the row structures of
	 *    those columns within the relaxed snode.
	 *    Note: The relaxed snodes are leaves of the supernodal etree,
	 *    therefore, the portion outside the rectangular supernode must be zero.
	 *
	 * Return value
	 * ============
	 *     0   success;
	 *    >0   number of bytes allocated when run out of memory.
	 *
	 */
	    GlobalLU_t Glu = pxgstrf_shared.Glu;
	    int i, k, ifrom, nextl, nsuper;
	    int          ito[] = new int[1];
	    int          krow, kmark, mem_error;
	    int[]        supno, lsub, xlsub, xlsub_end;

	    supno                 = Glu.supno;
	    xlsub                 = Glu.xlsub;
	    xlsub_end             = Glu.xlsub_end;
	    int[] nsuper_ = new int[1];
	    nsuper = NewNsuper(pnum, pxgstrf_shared, nsuper_);
	    Glu.nsuper = nsuper_[0];
	    Glu.xsup[nsuper]     = jcol;
	    Glu.xsup_end[nsuper] = kcol + 1;

	    nextl = 0;
	    for (i = jcol; i <= kcol; i++) {
		/* for each nonzero in A[*,i] */
		for (k = xa_begin[i]; k < xa_end[i]; k++) {
		    krow = asub[k];
		    kmark = marker[krow];
		    if ( kmark != kcol ) { /* First time visit krow */
			marker[krow] = kcol;
			col_lsub[nextl++] = krow;
		    }
	    	}
		supno[i] = nsuper;
	    }

	    if ( (mem_error = Glu_alloc(pnum, jcol, 2*nextl, LSUB, ito,
					pxgstrf_shared)) != 0 )
		return mem_error;

	    xlsub[jcol] = ito[0];
	    lsub        = Glu.lsub;
	    for (ifrom = 0; ifrom < nextl; ++ifrom)
		lsub[ito[0]++] = col_lsub[ifrom];
	    xlsub_end[jcol] = ito[0];

	    return 0;
	}

}

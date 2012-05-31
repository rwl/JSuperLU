package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_pdsp_defs.pxgstrf_shared_t;
import gov.lbl.superlu.Dlu_slu_mt_util.yes_no_t;
import gov.lbl.superlu.Dlu_supermatrix.NCPformat;
import gov.lbl.superlu.Dlu_supermatrix.SuperMatrix;

import static gov.lbl.superlu.Dlu.DEBUGlevel;
import static gov.lbl.superlu.Dlu_pdgstrf_pivotL.pdgstrf_pivotL;
import static gov.lbl.superlu.Dlu_pdgstrf_snode_bmod.pdgstrf_snode_bmod;
import static gov.lbl.superlu.Dlu_pdgstrf_snode_dfs.pdgstrf_snode_dfs;
import static gov.lbl.superlu.Dlu_pdutil.dprint_lu_col;
import static gov.lbl.superlu.Dlu_pmemory.Glu_alloc;
import static gov.lbl.superlu.Dlu_slu_mt_util.HICOL;
import static gov.lbl.superlu.Dlu_slu_mt_util.LOCOL;
import static gov.lbl.superlu.Dlu_slu_mt_util.MemType.LUSUP;





public class Dlu_pdgstrf_factor_snode {

	static
	int
	pdgstrf_factor_snode(
			     final int pnum,  /* process number */
			     final int jcol,
			     SuperMatrix A,
			     final double diag_pivot_thresh,
			     yes_no_t usepr[],
			     int    perm_r[],
			     int    inv_perm_r[], /* modified */
			     int    inv_perm_c[], /* in - used to find diagonal of Pc*A*Pc' */
			     int    xprune[],
			     int    marker[],
			     int    col_lsub[], /* values are irrevelant on entry
						  and on return */
			     double dense[],
			     double tempv[],
			     pxgstrf_shared_t pxgstrf_shared,
			     int    info[]
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
	 *
	 *   Factorize the artificial supernodes grouped at the bottom
	 *   of the etree.
	 *
	 */
	    GlobalLU_t   Glu = pxgstrf_shared.Glu;
	    int          singular;
	    NCPformat    Astore;
	    int kcol, icol, k, jsupno, fsupc, nsupr;
	    int ifrom, ito;
	    int          nextu, nextlu[] = new int[1];
	    int          pivrow[] = new int[1];
	    double       a[];
	    int[]        asub, xa_begin, xa_end, xusub, xusub_end,
	                 xsup, supno, xlusup, lsub, xlsub, xlsub_end;

	    lsub      = Glu.lsub;
	    xlsub     = Glu.xlsub;
	    xlsub_end = Glu.xlsub_end;
	    xusub     = Glu.xusub;
	    xusub_end = Glu.xusub_end;
	    xsup      = Glu.xsup;
	    supno     = Glu.supno;
	    xlusup    = Glu.xlusup;

	    singular = 0;
	    Astore   = (NCPformat) A.Store;
	    a        = Astore.nzval;
	    asub     = Astore.rowind;
	    xa_begin = Astore.colbeg;
	    xa_end   = Astore.colend;

	    kcol = jcol + pxgstrf_shared.pan_status[jcol].size;

	    /* Determine the union of the row structure of the supernode */
	    if ( (info[0] = pdgstrf_snode_dfs(pnum, jcol, kcol-1, asub, xa_begin, xa_end,
					   xprune, marker, col_lsub, pxgstrf_shared)) != 0 )
		return 0;

	    /*
	     * Factorize the relaxed supernode (jcol:kcol-1)
	     */
	    nextu        = Glu.nextu; /* xiaoye - race condition (no problem!) */
	    jsupno       = supno[jcol];
	    fsupc        = xsup[jsupno];
	    nsupr        = xlsub_end[fsupc] - xlsub[fsupc];
	    if ( (info[0] = Glu_alloc(pnum, jcol, nsupr*(kcol-jcol), LUSUP, nextlu,
				  pxgstrf_shared)) != 0 )
		return 0;

	    for (icol = jcol; icol < kcol; icol++) {
		xusub[icol] = xusub_end[icol] = nextu;
		xlusup[icol] = nextlu[0];

		/* Scatter into SPA dense[*] */
		for (k = xa_begin[icol]; k < xa_end[icol]; k++)
		    dense[asub[k]] = a[k];

		/* Numeric update within the supernode */
		pdgstrf_snode_bmod(pnum, icol, jsupno, fsupc, dense, tempv,
				   Glu, pxgstrf_shared.Gstat);

		if ( (info[0] = pdgstrf_pivotL
		                 (pnum, icol, diag_pivot_thresh, usepr, perm_r,
				  inv_perm_r, inv_perm_c, pivrow,
				  Glu, pxgstrf_shared.Gstat)) != 0 )
		    if ( singular == 0 ) singular = info[0];

		nextlu[0] += nsupr;

	if ( DEBUGlevel>= 2 ) {
	  if ( icol>=LOCOL && icol<=HICOL )
	    dprint_lu_col(pnum,"relax:",jcol,icol,kcol-jcol,pivrow[0],xprune,Glu);
	}

	    }

	    /* Store the row subscripts of kcol-1 for pruned graph */
	    k = ito = xlsub_end[jcol];
	    for (ifrom = xlsub[jcol]+kcol-jcol-1; ifrom < k; ++ifrom)
		lsub[ito++] = lsub[ifrom];
	    k = ito;
	    xprune[kcol-1] = k;
	    if (jcol < kcol-1) {    /* not a singleton */
		for (icol = jcol+1; icol < kcol; ++icol) xlsub_end[icol] = k;
		k = xlsub_end[jcol];
		xprune[jcol] = k;
		for (icol = jcol+1; icol < kcol; ++icol) xlsub[icol] = k;
	    }

	    info[0] = singular;
	    return 0;
	}

}

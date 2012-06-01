package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_pdsp_defs.pxgstrf_shared_t;

import static gov.lbl.superlu.Dlu_pxgstrf_synch.panel_t.RELAXED_SNODE;
import static gov.lbl.superlu.Dlu_await.await;

import static gov.lbl.superlu.Dlu_slu_mt_util.SUPER_FSUPC;
import static gov.lbl.superlu.Dlu_slu_mt_util.LOCOL;
import static gov.lbl.superlu.Dlu_slu_mt_util.HICOL;

import static gov.lbl.superlu.Dlu.DEBUGlevel;
import static gov.lbl.superlu.Dlu.printf;


public class Dlu_pxgstrf_mark_busy_descends {

	@SuppressWarnings("unused")
	static
	void
	pxgstrf_mark_busy_descends(int pnum, int jcol, int etree[],
				   pxgstrf_shared_t pxgstrf_shared,
				   int bcol[], int lbusy[])
	{
	/*
	 * -- SuperLU MT routine (version 1.0) --
	 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
	 * and Lawrence Berkeley National Lab.
	 * August 15, 1997
	 *
	 * Purpose
	 * =======
	 *
	 *   Mark busy panels in local "lbusy" array, used for linear pipelining.
	 *
	 *   When jcol begins, its busy descendant panels (if any) are bcol and
	 *   all the e-tree ancestors of bcol between bcol and jcol. This routine
	 *   marks those columns in the array lbusy, which is local to this
	 *   processor, to preserve a snapshot regardless of what the other
	 *   processors are doing meanwhile.
	 *
	 * Arguments
	 * =========
	 *
	 * jcol    (input) int
	 *         Current panel, with leading column jcol.
	 *
	 * etree   (input) int*
	 *         Elimination tree parent pointers.
	 *
	 * bcol    (input/output) int*
	 *         Farthest busy descendant of jcol.
	 *         On entry, it is the first column of the farthest busy panel.
	 *         On exit, it may be adjusted to the first column of the
	 *                  farthest busy supernode.
	 *
	 * lbusy   (input/output) int*
	 *         Initially all -1, lbusy(r) = jcol means that r was busy
	 *         at the beginning of computing jcol.
	 *
	 */
	    GlobalLU_t Glu = pxgstrf_shared.Glu;
	    int w,  kcol, fsupc, bcol_reg;
	    int xsup[];

	    bcol_reg = bcol[0];
	    if ( bcol_reg < jcol ) {

		/* -----------------------------------------------------------
		   Instead of waiting for the completion of "bcol", we can
		   pessimistically assume supno[bcol] == supno[bcol-1],
		   hence always mark as busy the supernode containing "bcol-1".
		   ----------------------------------------------------------- */
		if (pxgstrf_shared.pan_status[bcol_reg].type == RELAXED_SNODE) {
	if (false) {
		    if ( pxgstrf_shared.pan_status[bcol_reg].size < 0 )
		  	fsupc = bcol_reg + pxgstrf_shared.pan_status[bcol_reg].size;
		    else fsupc = bcol_reg;
	}
		    fsupc = bcol_reg;
		    w = pxgstrf_shared.pan_status[fsupc].size;
		    bcol_reg += w;
		    for (kcol = fsupc; kcol < bcol_reg; ++kcol)
			lbusy[kcol] = jcol;
		} else {
		    /* Find leading column "fsupc" in the supernode that
		       contains column "bcol-1" */
	if (false) {
		    if ( pxgstrf_shared.spin_locks[bcol_reg] != 0 ) /* WORSE PERFORMANCE!! */
			await( pxgstrf_shared.spin_locks[bcol_reg] );
	}
		    xsup = Glu.xsup;
		    fsupc = SUPER_FSUPC( xsup, Glu.supno[bcol_reg-1] );
		    for (kcol = fsupc; kcol < bcol_reg; ++kcol)	lbusy[kcol] = jcol;
		}

	if ( DEBUGlevel>=1 ) {
	if (jcol >= LOCOL && jcol <= HICOL)
	    printf("(%d) mark_busy_descends[1] jcol %d, bcol_reg %d, fsupc %d\n",
	           pnum, jcol, bcol_reg, fsupc);
	}

		/* Mark as busy all columns on the path between bcol_reg and jcol */
		for (kcol = bcol_reg; kcol < jcol; kcol = etree[kcol]) {
		    lbusy[kcol] = jcol;
		}

		/* INVARIANT: *bcol must be the first column of the farthest
		   busy supernode */
		bcol[0] = fsupc;

	    } /* if bcol_reg < jcol */
	}

}

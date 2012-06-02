/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 */
package gov.lbl.superlu;

import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_ABORT;
import static gov.lbl.superlu.Dlu_slu_mt_util.EMPTY;
import static gov.lbl.superlu.Dlu_util.ifill;
import static gov.lbl.superlu.Dlu_superlu_timer.SuperLU_timer_;

import static gov.lbl.superlu.Dlu_pxgstrf_synch.lu_locks_t.ULOCK;
import static gov.lbl.superlu.Dlu_pxgstrf_synch.lu_locks_t.LLOCK;
import static gov.lbl.superlu.Dlu_pxgstrf_synch.lu_locks_t.LULOCK;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_pdsp_defs.pxgstrf_shared_t;
import gov.lbl.superlu.Dlu_slu_mt_util.Gstat_t;
import gov.lbl.superlu.Dlu_slu_mt_util.MemType;

import static gov.lbl.superlu.Dlu.PROFILE;
import static gov.lbl.superlu.Dlu.stderr;
import static gov.lbl.superlu.Dlu.fprintf;
import static gov.lbl.superlu.Dlu.exit;


public class Dlu_pmemory {

	public static void
	XPAND_HINT(String memtype, int new_next, int jcol, int param) {
		fprintf(stderr, "Storage for %12s exceeded; Current column %d; Need at least %d;\n",
		memtype, jcol, new_next);
		fprintf(stderr, "You may set it by the %d-th parameter in routine sp_ienv().\n", param);
		SUPERLU_ABORT("Memory allocation failed");
	}

	/*
	 * Set up pointers for integer working arrays.
	 */
	public static void
	pxgstrf_SetIWork(int n, int panel_size, /*int iworkptr[], */int segrep[][],
			 int parent[][], int xplore[][], int repfnz[][], int panel_lsub[][],
			 int marker[][], int lbusy[][])
	{
	    segrep[0] = new int[n];                   /* n  */
	    parent[0] = new int[n];                   /* n  */
	    xplore[0] = new int[2*n];                 /* 2*n */
	    repfnz[0] = new int[panel_size*n];        /* w*n */
	    panel_lsub[0] = new int[panel_size*n];    /* w*n */
	    marker[0] = new int[3*n];                 /* 3*n */
	    lbusy[0]  = new int[n];                   /* n   */
	    ifill (repfnz[0], n * panel_size, EMPTY);

//	    segrep[0] = iworkptr;                                      /* n  */
//	    parent[0] = iworkptr + n;                                  /* n  */
//	    xplore[0] = iworkptr + 2*n;                                /* 2*n */
//	    repfnz[0] = iworkptr + 4*n;                                /* w*n */
//	    panel_lsub[0] = iworkptr + 4*n + panel_size*n;             /* w*n */
//	    marker[0] = iworkptr + 4*n + 2*panel_size*n;               /* 3*n */
//	    lbusy[0]  = iworkptr + (4+NO_MARKER)*n + 2*panel_size*n;   /* n   */
//	    ifill (repfnz[0], n * panel_size, EMPTY);
	}


	public static
	void
	copy_mem_int(int howmany, int old[], int new_[])
	{
	    int i;
	    int iold[] = old;
	    int inew[] = new_;
	    for (i = 0; i < howmany; i++) inew[i] = iold[i];
	}


	public static
	void
	user_bcopy(char src[], char dest[], int bytes)
	{
		System.arraycopy(src, 0, dest, 0, bytes);
//	    char s_ptr[], d_ptr[];
//
//	    s_ptr = src + bytes - 1;
//	    d_ptr = dest + bytes - 1;
//	    for (; d_ptr >= dest; --s_ptr, --d_ptr ) *d_ptr = *s_ptr;
	}



	@SuppressWarnings("unused")
	public static
	int[] intMalloc(int n)
	{
	    int buf[];
	    if ( (buf = new int [n]) == null ) {
		fprintf(stderr, "SUPERLU_MALLOC failed for buf in intMalloc()\n");
		exit (1);
	    }
	    return (buf);
	}

	@SuppressWarnings("unused")
	public static
	int[] intCalloc(int n)
	{
	    int buf[];
	    int i;
	    if ( (buf = new int [n]) == null ) {
		fprintf(stderr, "SUPERLU_MALLOC failed for buf in intCalloc()\n");
		exit (1);
	    }
	    for (i = 0; i < n; ++i) buf[i] = 0;
	    return (buf);
	}

	/*
	 * Allocate n elements storage from a global array.
	 * It uses lock for mutually exclusive access to the next position, so that
	 * more than one processors can call aalloc on the same array correctly.
	 * Return value: 0 - success
	 *              >0 - number of bytes allocated when run out of space
	 */
	public static
	int
	Glu_alloc(
		  final int pnum,     /* process number */
		  final int jcol,
		  final int num,      /* number of elements requested */
		  final MemType mem_type,
	          int   prev_next[],   /* return "next" value before allocation */
		  pxgstrf_shared_t pxgstrf_shared
		  )
	{
	    GlobalLU_t Glu = pxgstrf_shared.Glu;
	    Gstat_t    Gstat = pxgstrf_shared.Gstat;
	    int fsupc, nextl, nextu, new_next;
	    double   t = 0;

	    switch ( mem_type ) {

	      case LUSUP:
		/* Storage for the H-supernode is already set aside, so we do
		   not use lock here. */
		if ( Glu.map_in_sup[jcol] < 0 )
		    fsupc = jcol + Glu.map_in_sup[jcol];
		else fsupc = jcol;
		prev_next[0] = Glu.map_in_sup[fsupc];
		Glu.map_in_sup[fsupc] += num;

//	if (false) {
//		{
//		    int i, j;
//		    i = fsupc + part_super_h[fsupc];
//		    if ( Glu.dynamic_snode_bound == YES.ordinal() )
//		      new_next = Glu.nextlu;
//		    else new_next = Glu.map_in_sup[i];
//		    if (new_next < 0) { /* relaxed s-node */
//			for (i = fsupc+1; Glu.map_in_sup[i] < 0; ++i) ;
//			new_next = Glu.map_in_sup[i];
//		    }
//		    if ( Glu.map_in_sup[fsupc]>Glu.nzlumax
//			|| Glu.map_in_sup[fsupc]>new_next ) {
//			printf("(%d) jcol %d, map_[%d]=%d, map_[%d]=new_next %d\n",
//			       pnum, jcol, fsupc, Glu.map_in_sup[fsupc],
//			       i, new_next);
//			printf("(%d) snode type %d,size %d, |H-snode| %d\n",
//			       pnum, Glu.pan_status[fsupc].type,
//			       Glu.pan_status[fsupc].size, part_super_h[fsupc]);
//			for (j = fsupc; j < i; j += part_super_h[j])
//			    printf("(%d) H snode %d, size %d\n",
//				   pnum, j, part_super_h[j]);
//			SUPERLU_ABORT("LUSUP exceeded.");  /* xiaoye */
//		    }
//		}
//	}
		break;


	      case UCOL: case USUB:

	if (PROFILE) {
		t = SuperLU_timer_();
	}

		synchronized ( pxgstrf_shared.lu_locks[ULOCK.ordinal()] )
		{
		    nextu = Glu.nextu;
		    new_next = nextu + num;
		    if ( new_next > Glu.nzumax ) {
		        XPAND_HINT("U columns", new_next, jcol, 7);
		    }
		    prev_next[0] = nextu;
		    Glu.nextu = new_next;

		} /* end of critical region */

	if (PROFILE) {
		Gstat.procstat[pnum].cs_time += SuperLU_timer_() - t;
	}

		break;


		case LSUB:

	if (PROFILE) {
		t = SuperLU_timer_();
	}

	    synchronized( pxgstrf_shared.lu_locks[LLOCK.ordinal()] )
		{
		  nextl = Glu.nextl;
		  new_next = nextl + num;
		  if ( new_next > Glu.nzlmax ) {
		      XPAND_HINT("L subscripts", new_next, jcol, 8);
		  }
		  prev_next[0] = nextl;
		  Glu.nextl = new_next;

		} /* end of #pragama critical lock() */

	if (PROFILE) {
		Gstat.procstat[pnum].cs_time += SuperLU_timer_() - t;
	}

		  break;

	    }

	    return 0;
	}

	/*
	 * Dynamically set up storage image in lusup[*], using the supernode
	 * boundaries in H.
	 */
	public static
	int
	DynamicSetMap(
		      final int pnum,      /* process number */
		      final int jcol,      /* leading column of the s-node in H */
		      final int num,       /* number of elements requested */
		      pxgstrf_shared_t pxgstrf_shared
		      )
	{
		double t = 0;
	    GlobalLU_t Glu = pxgstrf_shared.Glu;
	    Gstat_t    Gstat = pxgstrf_shared.Gstat;
	    int nextlu, new_next;
	    int map_in_sup[] = Glu.map_in_sup; /* modified; memory mapping function */

	if (PROFILE) {
	    t = SuperLU_timer_();
	}

		synchronized ( pxgstrf_shared.lu_locks[LULOCK.ordinal()] )
	    {
		nextlu = Glu.nextlu;
		map_in_sup[jcol] = nextlu;
		new_next = nextlu + num;
		if ( new_next > Glu.nzlumax ) {
		    XPAND_HINT("L supernodes", new_next, jcol, 6);
		}
		Glu.nextlu = new_next;
	    } /* end of critical region */

	if (PROFILE) {
	    Gstat.procstat[pnum].cs_time += SuperLU_timer_() - t;
	}

	    return 0;
	}


}

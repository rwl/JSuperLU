/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 */
package lu.jsuper;

public class Dlu_pxgstrf_synch {

	/* Structure for the globally shared work queue */

	static class queue_t {
	    int       head, tail, count;
	    int       queue[];
	}

	enum lu_locks_t {
	    ULOCK,       /* locked once per column */
	    LLOCK,       /* locked once per supernode */
	    LULOCK,      /* locked once per column in L-supernode */
	    NSUPER_LOCK, /* locked once per supernode */
	    SCHED_LOCK,  /* locked once per panel, if succeeded each time */
	    NO_GLU_LOCKS
	}

	enum panel_t {
	    RELAXED_SNODE,
	    TREE_DOMAIN,   /* domain */
	    REGULAR_PANEL  /* non-domain */
	}

	enum pipe_state_t {
	    DONE,
	    BUSY,
	    CANGO,
	    CANPIPE,
	    UNREADY
	}

	static class pan_status_t {
	    panel_t      type;  /* panel type: 0 -- relaxed, also domain
				               1 -- domain
				               2 -- regular, non-domain */
	    pipe_state_t state; /* one of the 5 states in which the panel can be */
	    int          size;  /* in the leading column, the panel size is stored;
		                   in the other columns, the offset (negative)
			           to the leading column is stored */
	    int          ukids; /* number of kids not yet finished
				 * In linear pipeline --
				 *   if ukids[firstcol] = 0 then
				 *      the panel becomes a leaf (CANGO)
				 *   if ukids[firstcol] = 1 then
				 *      the panel can be taken as CANPIPE
				 */
	}


	/* The structure to record a relaxed supernode. */
	static class pxgstrf_relax_t {
	    int fcol;    /* first column of the relaxed supernode */
	    int size;    /* size of the relaxed supernode */
	}

}

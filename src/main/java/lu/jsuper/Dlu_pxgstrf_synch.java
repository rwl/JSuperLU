/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 */
package lu.jsuper;

import lu.jsuper.Dlu_pdsp_defs.pxgstrf_shared_t;
import lu.jsuper.Dlu_slu_mt_util.Gstat_t;
import lu.jsuper.Dlu_slu_mt_util.superlumt_options_t;
import static lu.jsuper.Dlu_slu_mt_util.SUPERLU_MAX;
import static lu.jsuper.Dlu_slu_mt_util.SUPERLU_MIN;
import static lu.jsuper.Dlu_slu_mt_util.SUPERLU_ABORT;
import static lu.jsuper.Dlu_slu_mt_util.EMPTY;
import static lu.jsuper.Dlu_slu_mt_util.STATE;
import static lu.jsuper.Dlu_slu_mt_util.Branch;

import static lu.jsuper.Dlu.printf;
import static lu.jsuper.Dlu.PRNTlevel;
import static lu.jsuper.Dlu.fprintf;
import static lu.jsuper.Dlu.stderr;
import static lu.jsuper.Dlu.PROFILE;
import static lu.jsuper.Dlu.DOMAINS;
import static lu.jsuper.Dlu.stdout;
import static lu.jsuper.Dlu.fflush;
import static lu.jsuper.Dlu.PREDICT_OPT;

import static lu.jsuper.Dlu_superlu_timer.SuperLU_timer_;


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

	static boolean SPLIT_TOP = false;

	static
	int
	ParallelInit(int n, pxgstrf_relax_t pxgstrf_relax[],
		     superlumt_options_t superlumt_options,
		     pxgstrf_shared_t pxgstrf_shared)
	{
	    int      etree[] = superlumt_options.etree;
	    int w, dad, ukids, i, j, k, rs, panel_size, relax;
	    int P, w_top, do_split = 0;
	    panel_t panel_type;
	    int      panel_histo[] = pxgstrf_shared.Gstat.panel_histo;
	    int nthr, concurrency, info;
	    Gstat_t Gstat = pxgstrf_shared.Gstat;

	    pxgstrf_shared.lu_locks = new pthread_mutex_t[lu_locks_t.NO_GLU_LOCKS.ordinal()];
	    for (i = 0; i < lu_locks_t.NO_GLU_LOCKS.ordinal(); ++i)
		pthread_mutex_init(pxgstrf_shared.lu_locks[i], null);

	if ( PRNTlevel==1 ) {
	    printf(".. ParallelInit() ... nprocs %2d\n", superlumt_options.nprocs);
	}

	    pxgstrf_shared.spin_locks = intCalloc(n);
	    pxgstrf_shared.pan_status = new pan_status_t[n+1];
	    pxgstrf_shared.fb_cols    = intMalloc(n+1);

	    panel_size = superlumt_options.panel_size;
	    relax = superlumt_options.relax;
	    w = SUPERLU_MAX(panel_size, relax) + 1;
	    for (i = 0; i < w; ++i) panel_histo[i] = 0;
	    pxgstrf_shared.num_splits = 0;

	    if ( (info = queue_init(pxgstrf_shared.taskq, n)) ) {
		fprintf(stderr, "ParallelInit(): %d\n", info);
		SUPERLU_ABORT("queue_init fails.");
	    }

	    /* Count children of each node in the etree. */
	    for (i = 0; i <= n; ++i) pxgstrf_shared.pan_status[i].ukids = 0;
	    for (i = 0; i < n; ++i) {
		dad = etree[i];
		++pxgstrf_shared.pan_status[dad].ukids;
	    }


	    /* Find the panel partitions and initialize each panel's status */

	if (PROFILE) {
	    Gstat.num_panels = 0;
	}

	    pxgstrf_shared.tasks_remain = 0;
	    rs = 1;   /* index for the next relaxed s-node */
	    w_top = panel_size/2;
	    if ( w_top == 0 ) w_top = 1;
	    P = 12;

	    for (i = 0; i < n; ) {
		if ( pxgstrf_relax[rs].fcol == i ) {
		    w = pxgstrf_relax[rs++].size;
		    panel_type = panel_t.RELAXED_SNODE;
		    pxgstrf_shared.pan_status[i].state = pipe_state_t.CANGO;
		} else {
		    /* Adjust panel_size so that a panel won't overlap with
		       the next relaxed snode.     */
	if (false) {
		    /* Only works when etree is postordered. */
		    w = SUPERLU_MIN(panel_size, pxgstrf_relax[rs].fcol - i);
	} else {
		    w = panel_size;
		    for (k = i + 1; k < SUPERLU_MIN(i + panel_size, n); ++k)
			if ( k == pxgstrf_relax[rs].fcol ) {
			    w = k - i;  /* panel stops at column k-1 */
			    break;
			}
		    if ( k == n ) w = n - i;
	}

	if (SPLIT_TOP) {
		    if ( do_split == 0 ) {
		  	if ( (n-i) < panel_size * P ) do_split = 1;
		    }
		    if ( do_split != 0 && w > w_top ) { /* split large panel */
		    	w = w_top;
		    	++pxgstrf_shared.num_splits;
		    }
	}
		    for (j = i+1; j < i + w; ++j)
			/* Do not allow panel to cross a branch point in the etree. */
			if ( pxgstrf_shared.pan_status[j].ukids > 1 ) break;
		    w = j - i;    /* j should start a new panel */
		    panel_type = panel_t.REGULAR_PANEL;
		    pxgstrf_shared.pan_status[i].state = pipe_state_t.UNREADY;
	if (DOMAINS) {
		    if ( in_domain[i] == TREE_DOMAIN ) panel_type = TREE_DOMAIN;
	}
		}

		if ( panel_type == panel_t.REGULAR_PANEL ) {
		    ++pxgstrf_shared.tasks_remain;
		    /*printf("nondomain panel %6d -- %6d\n", i, i+w-1);
		    fflush(stdout);*/
		}

		ukids = k = 0;
		for (j = i; j < i + w; ++j) {
		    pxgstrf_shared.pan_status[j].size = k--;
		    pxgstrf_shared.pan_status[j].type = panel_type;
		    ukids += pxgstrf_shared.pan_status[j].ukids;
		}
		pxgstrf_shared.pan_status[i].size = w; /* leading column */
		/* only count those kids outside the panel */
		pxgstrf_shared.pan_status[i].ukids = ukids - (w-1);
		panel_histo[w]++;

	if (PROFILE) {
		Gstat.panstat[i].size = w;
		++Gstat.num_panels;
	}

		pxgstrf_shared.fb_cols[i] = i;
		i += w;    /* move to the next panel */

	    } /* for i ... */

	    /* Dummy root */
	    pxgstrf_shared.pan_status[n].size = 1;
	    pxgstrf_shared.pan_status[n].state = pipe_state_t.UNREADY;

	if ( PRNTlevel==1 ) {
	    printf(".. Split: P %d, #nondomain panels %d\n", P, pxgstrf_shared.tasks_remain);
	}
	if (DOMAINS) {
	    EnqueueDomains(&pxgstrf_shared.taskq, list_head, pxgstrf_shared);
	} else {
	    EnqueueRelaxSnode(&pxgstrf_shared.taskq, n, pxgstrf_relax, pxgstrf_shared);
	}
	if ( PRNTlevel==1 ) {
	    printf(".. # tasks %d\n", pxgstrf_shared.tasks_remain);
	    fflush(stdout);
	}

	if (PREDICT_OPT) {
	    /* Set up structure describing children */
	    for (i = 0; i <= n; cp_firstkid[i++] = EMPTY);
	    for (i = n-1; i >= 0; i--) {
		dad = etree[i];
		cp_nextkid[i] = cp_firstkid[dad];
		cp_firstkid[dad] = i;
	    }
	}

	    return 0;
	} /* ParallelInit */


	/*
	 * Free the storage used by the parallel scheduling algorithm.
	 */
	static
	int ParallelFinalize(pxgstrf_shared_t pxgstrf_shared)
	{
	    /* Destroy mutexes */
	    int i;
	    for (i = 0; i < lu_locks_t.NO_GLU_LOCKS.ordinal(); ++i)
	        pthread_mutex_destroy( pxgstrf_shared.lu_locks[i] );

	    pxgstrf_shared.lu_locks = null;
	    pxgstrf_shared.spin_locks = null;
	    pxgstrf_shared.pan_status = null;
	    pxgstrf_shared.fb_cols = null;
	    pxgstrf_shared.Glu.map_in_sup = null;
	    queue_destroy(pxgstrf_shared.taskq);

	if ( PRNTlevel==1 ) {
	    printf(".. # panel splittings %d\n", pxgstrf_shared.num_splits);
	}

	    return 0;
	}

	int queue_init(queue_t q, int n)
	{
	    if ( n < 1 ) return (-1);

	    q.queue = new qitem_t[n];
	    q.count = 0;
	    q.head = 0;
	    q.tail = 0;

	    return 0;
	}

	int queue_destroy(queue_t q)
	{
	    q.queue = null;
	    return 0;
	}

	/*
	 * Return value: number of items in the queue
	 */
	int Enqueue(queue_t q, qitem_t item)
	{
	    q.queue[q.tail++] = item;
	    ++q.count;
	    return (q.count);
	}

	/*
	 * Return value: >= 0 number of items in the queue
	 *               = -1 queue is empty
	 */
	int Dequeue(queue_t q, qitem_t[] item)
	{
	    if ( q.count <= 0 ) return EMPTY;

	    item[0] = q.queue[q.head++];
	    --q.count;
	    return (q.count);
	}

	int QueryQueue(queue_t q)
	{
	    int     i;
	    printf("Queue count: %d\n", q.count);
	    for (i = q.head; i < q.tail; ++i)
		printf("%8d\titem %8d\n", i, q.queue[i]);

	    return 0;
	}

	int EnqueueRelaxSnode(queue_t q, int n, pxgstrf_relax_t pxgstrf_relax[],
			      pxgstrf_shared_t pxgstrf_shared)
	{
	    int rs, j, m;

	    m = pxgstrf_relax[0].size;
	    for (rs = 1; rs <= m; ++rs) {
		j = pxgstrf_relax[rs].fcol;
		q.queue[q.tail++] = j;
		q.count++;
		++pxgstrf_shared.tasks_remain;
	    }
	if ( PRNTlevel==1 ) {
	    printf(".. EnqueueRelaxSnode(): count %d\n", q.count);
	}
	    return 0;
	}

	/*
	 * Enqueue the initial independent domains.
	 * A pair of two numbers {root, fst_desc} is added in the queue.
	 */
	/*int EnqueueDomains(int P, queue_t *q, struct Branch **proc_domains_h)*/
	int EnqueueDomains(queue_t q, Branch list_head,
			   pxgstrf_shared_t pxgstrf_shared)
	{
	    Branch b, thrash;

	/*    for (pnum = 0; pnum < P; ++pnum) {
		for (b = proc_domains_h[pnum]; b != NULL; ) {*/
	    b = list_head;
	    while ( b != null ) {
		thrash = b;
		q.queue[q.tail++] = b.root;
		q.queue[q.tail++] = b.first_desc;
		q.count = q.count + 2;
//		STATE ( b.root ) = CANGO;
		pxgstrf_shared.pan_status[b.root].state = pipe_state_t.CANGO;
		++pxgstrf_shared.tasks_remain;
		b = b.next;
	    }
	    printf("EnqueueDomains(): count %d\n", q.count);
	    return 0;
	}

	static
	int NewNsuper(final int pnum, pxgstrf_shared_t pxgstrf_shared, int data[])
	{
	    int i;
	    double t;
	    mutex_t lock[] = &pxgstrf_shared.lu_locks[lu_locks_t.NSUPER_LOCK.ordinal()];
	    Gstat_t Gstat = pxgstrf_shared.Gstat;

	if (PROFILE) {
	    t = SuperLU_timer_();
	}

	    pthread_mutex_lock(lock);
	    {
	      i = ++(data[0]);
	    }
	    pthread_mutex_unlock(lock);

	if (PROFILE) {
	    Gstat.procstat[pnum].cs_time += SuperLU_timer_() - t;
	}

	    return i;
	}

	int lockon(int block[])
	{
	    while ( block[0] != 0 ) ; /* spin-wait */
	    block[0] = 1;
	    return 0;
	}

	int lockoff(int block[])
	{
	    block[0] = 0;
	    return 0;
	}


	}

}
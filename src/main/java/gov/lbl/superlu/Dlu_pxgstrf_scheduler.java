package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.pxgstrf_shared_t;
import gov.lbl.superlu.Dlu_slu_mt_util.Gstat_t;

import static gov.lbl.superlu.Dlu_slu_mt_util.EMPTY;
import static gov.lbl.superlu.Dlu_slu_mt_util.STATE;
import static gov.lbl.superlu.Dlu_slu_mt_util.DADPANEL;
import static gov.lbl.superlu.Dlu_slu_mt_util.TIC;
import static gov.lbl.superlu.Dlu_slu_mt_util.TOC;
import static gov.lbl.superlu.Dlu_slu_mt_util.how_selected_t.DADPAN;
import static gov.lbl.superlu.Dlu_slu_mt_util.how_selected_t.PIPE;
import static gov.lbl.superlu.Dlu_slu_mt_util.how_selected_t.NOPIPE;
import static gov.lbl.superlu.Dlu_pxgstrf_synch.lu_locks_t.SCHED_LOCK;
import static gov.lbl.superlu.Dlu_pxgstrf_synch.pipe_state_t.BUSY;
import static gov.lbl.superlu.Dlu_pxgstrf_synch.pipe_state_t.CANGO;

import static gov.lbl.superlu.Dlu.DOMAINS;
import static gov.lbl.superlu.Dlu.PROFILE;
import static gov.lbl.superlu.Dlu.DEBUG;
import static gov.lbl.superlu.Dlu.POSTORDER;
import static gov.lbl.superlu.Dlu.printf;

import static gov.lbl.superlu.Dlu_superlu_timer.SuperLU_timer_;


public class Dlu_pxgstrf_scheduler {

	static
	void
	pxgstrf_scheduler(final int pnum, final int n, final int etree[],
			  int cur_pan[], int bcol[], pxgstrf_shared_t pxgstrf_shared)
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
	 * pxgstrf_scheduler() gets a panel for the processor to work on.
	 * It schedules a panel in decreasing order of priority:
	 *   (1) the current panel's parent, if it can be done without pipelining
	 *   (2) any other panel in the queue that can be done without pipelining
	 *       ("CANGO" status)
	 *   (3) any other panel in the queue that can be done with pipelining
	 *       ("CANPIPE" status)
	 *
	 * Arguments
	 * =========
	 * pnum    (input) int
	 *         Processor number.
	 *
	 * n       (input) int
	 *         Column dimension of the matrix.
	 *
	 * etree   (input) int*
	 *         Elimination tree of A'*A, size n.
	 *         Note: etree is a vector of parent pointers for a forest whose
	 *         vertices are the integers 0 to n-1; etree[root] = n.
	 *
	 * cur_pan (input/output) int*
	 *         On entry, the current panel just finished by this processor;
	 *         On exit, [0, n-1]: the new panel to work on;
	 *                  EMPTY:    failed to get any work, will try later;
	 *                  n:        all panels are taken; ready to terminate.
	 *
	 * taskq   (input/output) queue_t*
	 *         Global work queue.
	 *
	 * fb_cols (input/output) int*
	 *         The farthest busy descendant of each (leading column of the) panel.
	 *
	 * bcol    (output) int*
	 *         The most distant busy descendant of cur_pan in the *linear*
	 *         pipeline of busy descendants. If all proper descendants of
	 *         cur_pan are done, bcol is returned equal to cur_pan.
	 *
	 * Defining terms
	 * ==============
	 *   o parent(panel) = etree(last column in the panel)
	 *   o the kids of a panel = collective kids of all columns in the panel
	 *     kids[REP] = SUM_{j in panel} ( kids[j] )
	 *   o linear pipeline - what does it mean in the panel context?
	 *       if ukids[REP] = 0, then the panel becomes a leaf (CANGO)
	 *       if ukids[REP] = 1 && ukids[firstcol] = 1, then the panel can
	 *                       be taken with pipelining (CANPIPE)
	 *
	 * NOTES
	 * =====
	 *   o When a "busy" panel finishes, if its parent has only one remaining
	 *     undone child there is no check to see if the parent should change
	 *     from "unready" to "canpipe". Thus a few potential pipelinings will
	 *     be lost, but checking out this pipeline opportunity may be costly.
	 *
	 */

	    int dad, dad_ukids, jcol, w, j;
	    int fb_cols[] = pxgstrf_shared.fb_cols;
	    queue_t taskq[] = &pxgstrf_shared.taskq;
	    Gstat_t Gstat = pxgstrf_shared.Gstat;
	    double t;

	    jcol = cur_pan[0];
	    if ( jcol != EMPTY ) {
	if (DOMAINS) {
		if ( in_domain[jcol] == TREE_DOMAIN )
		    dad = etree[jcol];
		else
	}
		    dad = DADPANEL (jcol);
	    }

	    /* w_top = sp_ienv(1)/2;
	       if ( w_top == 0 ) w_top = 1;*/

	if (PROFILE) {
	    TIC(t);
	}
	    pthread_mutex_lock( pxgstrf_shared.lu_locks[SCHED_LOCK.ordinal()] );

	{   /* ---- START CRITICAL SECTION ---- */

	    /* Update the status of the current panel and its parent, so that
	     * the other processors waiting on it can proceed.
	     * If all siblings are done, and dad is not busy, then take dad.
	     */
	    if ( jcol != EMPTY ) { /* jcol was just finished by this processor */
		dad_ukids = --pxgstrf_shared.pan_status[dad].ukids;

	if (DEBUG) {
		printf("(%d) DONE %d in Scheduler(), dad %d, STATE %d, dad_ukids %d\n",
		       pnum, jcol, dad, STATE(dad), dad_ukids);
	}

		if ( dad_ukids == 0 && STATE( dad ) > BUSY.ordinal() ) { /* dad not started */
		    jcol = dad;
	if (DEBUG) {
		    printf("(%d) Scheduler[1] Got dad %d, STATE %d\n",
			   pnum, jcol, STATE(dad));
	}
	if (PROFILE) {
		    ++(Gstat.panhows[DADPAN.ordinal()]);
	}
		} else {
		    /* Try to get a panel from the task Q. */
		    while ( true ) {
			/*>>if ( (j = Dequeue(taskq, &item)) == EMPTY ) {*/
			if ( taskq.count <= 0 ) {
			    jcol = EMPTY;
			    break;
			} else {
			    jcol = taskq.queue[taskq.head++];
			    --taskq.count;
			    if ( STATE( jcol ) >= CANGO.ordinal() ) { /* CANGO or CANPIPE */
	if (DEBUG) {
				printf("(%d) Dequeue[1] Got %d, STATE %d, Qcount %d\n",
				       pnum, jcol, STATE(jcol), j);
	}
	if (PROFILE) {
				if (STATE( jcol ) == CANGO.ordinal()) ++(Gstat.panhows[NOPIPE.ordinal()]);
				else ++(Gstat.panhows[PIPE.ordinal()]);
	}
			        break;
			    }
			}
		    } /* while */
		}
	    } else {
		/*
		 * jcol was EMPTY; Try to get a panel from the task Q.
		 */
	    	while ( true ) {
	    	    /*>>if ( (j = Dequeue(taskq, &item)) == EMPTY ) {*/
		    if ( taskq.count <= 0 ) {
			jcol = EMPTY;
			break;
		    } else {
			jcol = taskq.queue[taskq.head++];
			--taskq.count;
			if ( STATE( jcol ) >= CANGO.ordinal() ) { /* CANGO or CANPIPE */
	if (DEBUG) {
			    printf("(%d) Dequeue[2] Got %d, STATE %d, Qcount %d\n",
				   pnum, jcol, STATE(jcol), j);
	}
	if (PROFILE) {
			    if (STATE( jcol ) == CANGO.ordinal()) ++(Gstat.panhows[NOPIPE.ordinal()]);
			    else ++(Gstat.panhows[PIPE.ordinal()]);
	}
			    break;
			}
		    }
		} /* while */
	    }

	    /*
	     * Update the status of the new panel "jcol" and its parent "dad".
	     */
	    if ( jcol != EMPTY ) {
		    --pxgstrf_shared.tasks_remain;
	if (DOMAINS) {
		if ( in_domain[jcol] == TREE_DOMAIN ) {
		    /* Dequeue the first descendant of this domain */
		    *bcol = taskq.queue[taskq.head++];
		    --taskq.count;
		} else
	}
		{
		    STATE( jcol ) = BUSY;
		    w = pxgstrf_shared.pan_status[jcol].size;

		    for (j = jcol; j < jcol+w; ++j) pxgstrf_shared.spin_locks[j] = 1;
		    dad = DADPANEL (jcol);
		    if ( dad < n && pxgstrf_shared.pan_status[dad].ukids == 1 ) {
			STATE( dad ) = CANPIPE;
			/*>> j = Enqueue(taskq, dad);*/
			taskq.queue[taskq.tail++] = dad;
			++taskq.count;
	if (DEBUG) {
			printf("(%d) Enqueue() %d's dad %d .CANPIPE, Qcount %d\n",
			       pnum, jcol, dad, j);
	}
		    }

	if (PROFILE) {
		    Gstat.procstat[pnum].panels++;
	}

		    /* Find the farthest busy descendant of the new panel
		       and its parent.*/
		    bcol[0] = fb_cols[jcol];
	if (DEBUG) {
		    printf("(%d) Scheduler[2] fb_cols[%d]=%d, STATE %d\n",
			   pnum, jcol, bcol[0], STATE( bcol[0] ));
	}
		    while ( STATE( bcol[0] ) == DONE ) bcol[0] = DADPANEL (*bcol);
		    fb_cols[dad] = bcol[0];

		} /* else regular_panel */

	    } /* if jcol != empty */

	    cur_pan[0] = jcol;

	if (DEBUG) {
	    printf("(%d) Exit C.S. tasks_remain %d, cur_pan %d\n",
		   pnum, pxgstrf_shared.tasks_remain, jcol);
	}

	} /* ---- END CRITICAL SECTION ---- */

	    pthread_mutex_unlock( pxgstrf_shared.lu_locks[SCHED_LOCK.ordinal()] );

	if (PROFILE) {
	    Gstat.procstat[pnum].cs_time += SuperLU_timer_() - t;
	}

	    return;
	}


	/* Fix the order of the panels to be taken. */
	static
	void
	Preorder(final int pnum, final int n, final int etree[], int cur_pan[],
	         queue_t taskq[], int fb_cols[], int bcol[],
		 pxgstrf_shared_t pxgstrf_shared)
	{
	    int w, dad, dad_ukids;

	POSTORDER = false;
	if (POSTORDER) {
	    if ( cur_pan[0] == EMPTY ) {
		cur_pan[0] = 0;
	    } else {
		w = pxgstrf_shared.pan_status[cur_pan[0]].size;
		cur_pan[0] += w;
	    }
	} else { /* Breadth-first bottom up */
	    if ( cur_pan[0] != EMPTY ) {
		dad = DADPANEL (cur_pan[0]);
		dad_ukids = --pxgstrf_shared.pan_status[dad].ukids;
		if ( dad_ukids == 0 ) {
		    taskq.queue[taskq.tail++] = dad;
		    ++taskq.count;
		}
	    }
	    cur_pan[0] = taskq.queue[taskq.head++];
	    --taskq.count;
	}
	    --pxgstrf_shared.tasks_remain;
	    bcol[0] = cur_pan[0];
	}

}

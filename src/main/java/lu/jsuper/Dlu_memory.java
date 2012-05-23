/*! @file Dlu_memory.java
 * \brief Precision-independent memory-related routines
 *
 * <pre>
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 * </pre>
 */
package lu.jsuper;

import static lu.jsuper.Dlu_slu_util.ABORT;
import static lu.jsuper.Dlu_slu_util.EMPTY;
import static lu.jsuper.Dlu_util.ifill;


/** Precision-independent memory-related routines.
    (Shared by [sdcz]memory.java) **/
public class Dlu_memory {

	/**! \brief Set up pointers for integer working arrays.
	 */
	public static void SetIWork(int m, int n, int panel_size, int iworkptr[],
			int segrep[][], int parent[][], int xplore[][], int repfnz[][],
			int panel_lsub[][], int xprune[][], int marker[][])
	{
	    segrep[0] = iworkptr;
	    parent[0] = iworkptr + m;
	    xplore[0] = parent[0] + m;
	    repfnz[0] = xplore[0] + m;
	    panel_lsub[0] = repfnz[0] + panel_size * m;
	    xprune[0] = panel_lsub[0] + panel_size * m;
	    marker[0] = xprune[0] + n;
	    ifill (repfnz[0], m * panel_size, EMPTY);
	    ifill (panel_lsub[0], m * panel_size, EMPTY);
	}

	public static void copy_mem_int(int howmany, int old[], int new_[])
	{
	    int i;
	    int iold[] = old;
	    int inew[] = new_;
	    for (i = 0; i < howmany; i++) inew[i] = iold[i];
	}

	public static int[] intMalloc(int n) {
		int[] buf = null;
		try {
			buf = new int[n];
		} catch (OutOfMemoryError e) {
			ABORT("SUPERLU_MALLOC fails for buf in intMalloc()", e);
		}
		return (buf);
	}

	public static int[] intCalloc(int n) {
		int buf[] = null;
		int i;
		try {
			buf = new int[n];
		} catch (OutOfMemoryError e) {
			ABORT("SUPERLU_MALLOC fails for buf in intCalloc()", e);
		}
	    for (i = 0; i < n; ++i) buf[i] = 0;
	    return (buf);
	}

}

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


/** Precision-independent memory-related routines.
    (Shared by [sdcz]memory.java) **/
public class Dlu_memory {

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

/*! @file Dlu_superlu_timer.java
 * \brief Returns the time used
 *
 * <pre>
 * Purpose
 * =======
 *
 * Returns the time in seconds used by the process.
 *
 * Note: the timer function call is machine dependent. Use conditional
 *       compilation to choose the appropriate function.
 * </pre>
 */

package lu.jsuper;


public class Dlu_superlu_timer {

	public static double SuperLU_timer_() {
		return (double) System.nanoTime();
	}

}

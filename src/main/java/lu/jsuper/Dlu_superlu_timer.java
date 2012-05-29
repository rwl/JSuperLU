/*
 * -- SuperLU MT routine (alpha version) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * Purpose
 * =======
 *	Returns the time in seconds used by the process.
 *
 * Note: the timer function call is machine dependent. Use conditional
 *       compilation to choose the appropriate function.
 *
 */
package lu.jsuper;

public class Dlu_superlu_timer {

	static double SuperLU_timer_()
	{
	    double dclock = System.nanoTime();
	    return (dclock);
	}

	static
	double usertimer_()
	{
	    return SuperLU_timer_();
	}

}

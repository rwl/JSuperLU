/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * This routine should NOT be optimized.
 */
package lu.jsuper;

public class Dlu_await {

	@SuppressWarnings("unused")
	static
	int await(int status[])
	{
	    int i, j, k, randnum;

	    /* randnum = ( random() & 0xff ); */
	    randnum = 0;
	    while ( status[0] != 0 ) ;
	if (false) {
	    {
		/* Length better be adaptive to the number of processors */
		k = randnum;
		for (i = 0; i < randnum; ++i) {
		    j += k;
		    k = -k;
		}
	    }
	}
	    return 0;
	}

}

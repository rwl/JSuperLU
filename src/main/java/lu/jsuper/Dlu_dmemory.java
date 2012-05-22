package lu.jsuper;

import static lu.jsuper.Dlu_slu_util.ABORT;


public class Dlu_dmemory {

	public static double[] doubleMalloc(int n) {
		double[] buf = null;
		try {
			buf = new double[n];
		} catch (OutOfMemoryError e) {
			ABORT("SUPERLU_MALLOC failed for buf in doubleMalloc()\n", e);
		}
		return (buf);
	}

}

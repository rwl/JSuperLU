package lu.jsuper;

import static lu.jsuper.Dlu.printf;

import java.io.PrintStream;

public class Dlu {

	public static final PrintStream stdout = System.out;
	public static final PrintStream stderr = System.err;

	public static int PRNTlevel = 0;
	public static int DEBUGlevel = 0;
	public static boolean PROFILE = false;
	public static boolean PREDICT_OPT = false;
	public static boolean USE_VENDOR_BLAS = false;
	public static boolean GEMV2 = false;
	public static boolean SCATTER_FOUND = false;

	public static void printf(String format, Object... args) {
		System.out.printf(format, args);
	}

	public static void fprintf(PrintStream stream, String format, Object... args) {
		stream.printf(format, args);
	}

	public static void fflush(PrintStream stream) {
		stream.flush();
	}

	public static double fabs(double a) {
		return Math.abs(a);
	}

	public static void exit(int status) {
		System.exit(status);
	}

}

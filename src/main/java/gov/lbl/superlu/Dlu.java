package gov.lbl.superlu;

import java.io.PrintStream;

public class Dlu {

	public static final PrintStream stdout = System.out;
	public static final PrintStream stderr = System.err;

	public static boolean DEBUG = false;

	public static int PRNTlevel = 0;
	public static int DEBUGlevel = 0;
	public static boolean PROFILE = false;
	public static boolean PREDICT_OPT = false;
	public static boolean USE_VENDOR_BLAS = false;
	public static boolean GEMV2 = false;
	public static boolean SCATTER_FOUND = false;

	public static boolean CHK_COLORDER = false;
	public static boolean ZFD_PERM = false;
	public static boolean CHK_NZCNT = false;
	public static boolean DOMAINS = false;
	public static boolean CHK_EXPAND = false;
	public static boolean POSTORDER = false;
	public static boolean CHK_PIVOT = false;
	public static boolean CHK_DFS = false;
	public static boolean DOPRINT = false;
	public static boolean CHK_PRUNE = false;
	public static boolean COMPRESS_LUSUP = false;

	public static void printf(String format, Object... args) {
		System.out.printf(format, args);
	}

	public static void fprintf(PrintStream stream, String format, Object... args) {
		stream.printf(format, args);
	}

	public static void fflush(PrintStream stream) {
		stream.flush();
	}

	public static void fclose(PrintStream stream) {
		stream.close();
	}

	public static double fabs(double a) {
		return Math.abs(a);
	}

	public static void exit(int status) {
		System.exit(status);
	}

}

package gov.lbl.superlu;

import java.io.PrintStream;

import org.netlib.blas.BLAS;

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

	public static String getenv(String name) {
		return System.getenv(name);
	}

	public static int[] d2i(double[] d) {
		int[] iarray = new int[d.length];
		for (int i = 0; i < d.length; i++)
			iarray[i] = (int) d[i];
		return iarray;
	}

	public static void dtrsm(String side, String uplo, String transa, String diag,
			int m, int n, double alpha, double[] a, int a_offset, int lda,
			double[] b, int b_offset, int ldb) {

		BLAS blas = BLAS.getInstance();

		int k;
		double[] A, B;

		k = side.equalsIgnoreCase("L") ? m : n;

		A = new double[lda*k];
		System.arraycopy(a, a_offset, A, 0, lda*k);

		B = new double[ldb*n];
		System.arraycopy(b, b_offset, B, 0, ldb*n);

 		blas.dtrsm(side, uplo, transa, diag, m, n, alpha,
		       A, lda, B, ldb);

 		System.arraycopy(B, 0, b, b_offset, ldb*n);

	}

	public static void dgemm(String transa, String transb, int m, int n,
			int k, double alpha, double[] a, int a_offset, int lda,
			double[] b, int b_offset, int ldb, double beta, double[] c,
			int Ldc) {

		BLAS blas = BLAS.getInstance();

		int ka, kb;
		double[] A, B;

		ka = transa.equalsIgnoreCase("N") ? k : m;
		kb = transb.equalsIgnoreCase("N") ? n : k;

 		A = new double[lda*ka];
 		System.arraycopy(a, a_offset, A, 0, lda*ka);

 		B = new double[ldb*kb];
 		System.arraycopy(b, b_offset, B, 0, ldb*kb);

		blas.dgemm( "N", "N", m, n, k, alpha,
			A, lda, B, ldb, beta, c, Ldc );

	}

	public static void dtrsv(String uplo, String trans, String diag, int n,
			double[] a, int lda, double[] x, int incx) {
		dtrsv(uplo, trans, diag, n, a, 0, lda, x, 0, incx);
	}

	public static void dtrsv(String uplo, String trans, String diag, int n,
			double[] a, int a_offset, int lda, double[] x, int x_offset, int incx) {

		BLAS blas = BLAS.getInstance();

		double[] A, X;

		A = new double[lda*n];
		System.arraycopy(a, a_offset, A, 0, lda*n);

		X = new double[];
		System.arraycopy(x, x_offset, X, 0, length);

		blas.dtrsv(uplo, trans, diag, n, A, lda, X, incx);

		System.arraycopy(X, 0, x, x_offset, length);
	}

	public static void dgemv(String trans, int m, int n, double alpha,
			double[] a, int a_offset, int lda, double[] x, int x_offset, int incx,
			double beta, double[] y, int y_offset, int incy) {

		BLAS blas = BLAS.getInstance();

		double[] A, X, Y;

		blas.dgemv(trans, m, n, alpha, A, lda, X, incx, beta, Y, incy);
	}



}

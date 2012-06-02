/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 */
package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_slu_mt_util.Gstat_t;
import gov.lbl.superlu.Dlu_slu_mt_util.superlu_memusage_t;
import gov.lbl.superlu.Dlu_slu_mt_util.trans_t;
import gov.lbl.superlu.Dlu_supermatrix.DNformat;
import gov.lbl.superlu.Dlu_supermatrix.Dtype_t;
import gov.lbl.superlu.Dlu_supermatrix.Mtype_t;
import gov.lbl.superlu.Dlu_supermatrix.NCPformat;
import gov.lbl.superlu.Dlu_supermatrix.NCformat;
import gov.lbl.superlu.Dlu_supermatrix.SCPformat;
import gov.lbl.superlu.Dlu_supermatrix.SCformat;
import gov.lbl.superlu.Dlu_supermatrix.Stype_t;
import gov.lbl.superlu.Dlu_supermatrix.SuperMatrix;

import static gov.lbl.superlu.Dlu_slu_mt_util.FALSE;
import static gov.lbl.superlu.Dlu_slu_mt_util.TRUE;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_ABORT;
import static gov.lbl.superlu.Dlu_slu_mt_util.trans_t.NOTRANS;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_MAX;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.FACT;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.RELAX;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.SOLVE;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.ETREE;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.EQUIL;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.RCOND;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.REFINE;

import static gov.lbl.superlu.Dlu.printf;
import static gov.lbl.superlu.Dlu.fflush;
import static gov.lbl.superlu.Dlu.stdout;
import static gov.lbl.superlu.Dlu.fabs;

import static gov.lbl.superlu.Dlu_pdmemory.doubleMalloc;
import static gov.lbl.superlu.Dlu_pmemory.intCalloc;
import static gov.lbl.superlu.Dlu_pmemory.intMalloc;

import static gov.lbl.superlu.Dlu_dsp_blas3.sp_dgemm;


public class Dlu_pdutil {

	public static
	void
	dCreate_CompCol_Matrix(SuperMatrix A, int m, int n, int nnz, double nzval[],
			      int rowind[], int colptr[],
			      Stype_t stype, Dtype_t dtype, Mtype_t mtype)
	{
	    NCformat Astore;

	    A.Stype = stype;
	    A.Dtype = dtype;
	    A.Mtype = mtype;
	    A.nrow = m;
	    A.ncol = n;
	    A.Store = new NCformat();
	    Astore = (NCformat) A.Store;
	    Astore.nnz = nnz;
	    Astore.nzval = nzval;
	    Astore.rowind = rowind;
	    Astore.colptr = colptr;
	}

	public static
	void
	dCreate_CompCol_Permuted(SuperMatrix A, int m, int n, int nnz, double nzval[],
				 int rowind[], int colbeg[], int colend[],
				 Stype_t stype, Dtype_t dtype, Mtype_t mtype)
	{
	    NCPformat Astore;

	    A.Stype = stype;
	    A.Dtype = dtype;
	    A.Mtype = mtype;
	    A.nrow = m;
	    A.ncol = n;
	    A.Store = new NCPformat();
	    Astore = (NCPformat) A.Store;
	    Astore.nnz = nnz;
	    Astore.nzval = nzval;
	    Astore.rowind = rowind;
	    Astore.colbeg = colbeg;
	    Astore.colend = colend;
	}
	/*
	 * Convert a row compressed storage into a column compressed storage.
	 */
	public static
	void
	dCompRow_to_CompCol(int m, int n, int nnz,
			    double a[], int colind[], int rowptr[],
			    double at[][], int rowind[][], int colptr[][])
	{
	    int i, j, col, relpos;
	    int marker[];

	    /* Allocate storage for another copy of the matrix. */
	    at[0] = (double []) doubleMalloc(nnz);
	    rowind[0] = (int []) intMalloc(nnz);
	    colptr[0] = (int []) intMalloc(n+1);
	    marker = (int []) intCalloc(n);

	    /* Get counts of each column of A, and set up column pointers */
	    for (i = 0; i < m; ++i)
		for (j = rowptr[i]; j < rowptr[i+1]; ++j) ++marker[colind[j]];
	    (colptr[0])[0] = 0;
	    for (j = 0; j < n; ++j) {
		(colptr[0])[j+1] = (colptr[0])[j] + marker[j];
		marker[j] = (colptr[0])[j];
	    }

	    /* Transfer the matrix into the compressed column storage. */
	    for (i = 0; i < m; ++i) {
		for (j = rowptr[i]; j < rowptr[i+1]; ++j) {
		    col = colind[j];
		    relpos = marker[col];
		    (rowind[0])[relpos] = i;
		    (at[0])[relpos] = a[j];
		    ++marker[col];
		}
	    }
	}


	/* Copy matrix A into matrix B. */
	public static
	void
	dCopy_CompCol_Matrix(SuperMatrix A, SuperMatrix B)
	{
	    NCformat Astore, Bstore;
	    int      ncol, nnz, i;

	    B.Stype = A.Stype;
	    B.Dtype = A.Dtype;
	    B.Mtype = A.Mtype;
	    B.nrow  = A.nrow;;
	    B.ncol  = ncol = A.ncol;
	    Astore   = (NCformat) A.Store;
	    Bstore   = (NCformat) B.Store;
	    Bstore.nnz = nnz = Astore.nnz;
	    for (i = 0; i < nnz; ++i)
		((double [])Bstore.nzval)[i] = ((double [])Astore.nzval)[i];
	    for (i = 0; i < nnz; ++i) Bstore.rowind[i] = Astore.rowind[i];
	    for (i = 0; i <= ncol; ++i) Bstore.colptr[i] = Astore.colptr[i];
	}


	public static
	int dPrint_CompCol_Matrix(SuperMatrix A)
	{
	    NCformat     Astore;
	    int          i;
	    double       dp[];

	    printf("\nCompCol matrix: ");
	    printf("Stype %d, Dtype %d, Mtype %d\n", A.Stype,A.Dtype,A.Mtype);
	    Astore = (NCformat) A.Store;
	    dp = (double []) Astore.nzval;
	    printf("nrow %d, ncol %d, nnz %d\n", A.nrow,A.ncol,Astore.nnz);
	    printf("\nnzval: ");
	    for (i = 0; i < Astore.nnz; ++i) printf("%f  ", dp[i]);
	    printf("\nrowind: ");
	    for (i = 0; i < Astore.nnz; ++i) printf("%d  ", Astore.rowind[i]);
	    printf("\ncolptr: ");
	    for (i = 0; i <= A.ncol; ++i) printf("%d  ", Astore.colptr[i]);
	    printf("\nend CompCol matrix.\n");

	    return 0;
	}

	public static
	int dPrint_Dense_Matrix(SuperMatrix A)
	{
	    DNformat     Astore;
	    int          i;
	    double       dp[];

	    printf("\nDense matrix: ");
	    printf("Stype %d, Dtype %d, Mtype %d\n", A.Stype,A.Dtype,A.Mtype);
	    Astore = (DNformat) A.Store;
	    dp = (double []) Astore.nzval;
	    printf("nrow %d, ncol %d, lda %d\n", A.nrow,A.ncol,Astore.lda);
	    printf("\nnzval: ");
	    for (i = 0; i < A.nrow; ++i) printf("%f  ", dp[i]);
	    printf("\nend Dense matrix.\n");

	    return 0;
	}

	public static
	void
	dCreate_Dense_Matrix(SuperMatrix X, int m, int n, double x[], int ldx,
			    Stype_t stype, Dtype_t dtype, Mtype_t mtype)
	{
	    DNformat    Xstore;

	    X.Stype = stype;
	    X.Dtype = dtype;
	    X.Mtype = mtype;
	    X.nrow = m;
	    X.ncol = n;
	    X.Store = new DNformat();
	    Xstore = (DNformat) X.Store;
	    Xstore.lda = ldx;
	    Xstore.nzval = (double []) x;
	}

	public static
	void
	dCopy_Dense_Matrix(int M, int N, double X[], int ldx, double Y[], int ldy)
	{
	/*
	 *
	 *  Purpose
	 *  =======
	 *
	 *  Copies a two-dimensional matrix X to another matrix Y.
	 */
	    int    i, j;

	    for (j = 0; j < N; ++j)
	        for (i = 0; i < M; ++i)
	            Y[i + j*ldy] = X[i + j*ldx];
	}

	public static
	void
	dCreate_SuperNode_Matrix(SuperMatrix L, int m, int n, int nnz, double nzval[],
				int nzval_colptr[], int rowind[], int rowind_colptr[],
				int col_to_sup[], int sup_to_col[],
				Stype_t stype, Dtype_t dtype, Mtype_t mtype)
	{
	    SCformat Lstore;

	    L.Stype = stype;
	    L.Dtype = dtype;
	    L.Mtype = mtype;
	    L.nrow = m;
	    L.ncol = n;
	    L.Store = new SCformat();
	    Lstore = (SCformat) L.Store;
	    Lstore.nnz = nnz;
	    Lstore.nsuper = col_to_sup[n];
	    Lstore.nzval = nzval;
	    Lstore.nzval_colptr = nzval_colptr;
	    Lstore.rowind = rowind;
	    Lstore.rowind_colptr = rowind_colptr;
	    Lstore.col_to_sup = col_to_sup;
	    Lstore.sup_to_col = sup_to_col;

	}

	public static
	void
	dCreate_SuperNode_Permuted(SuperMatrix L, int m, int n, int nnz,
				   double nzval[],
				   int nzval_colbeg[], int nzval_colend[],
				   int rowind[], int rowind_colbeg[], int rowind_colend[],
				   int col_to_sup[],
				   int sup_to_colbeg[], int sup_to_colend[],
				   Stype_t stype, Dtype_t dtype, Mtype_t mtype)
	{
	    SCPformat Lstore;

	    L.Stype = stype;
	    L.Dtype = dtype;
	    L.Mtype = mtype;
	    L.nrow = m;
	    L.ncol = n;
	    L.Store = new SCPformat();
	    Lstore = (SCPformat) L.Store;
	    Lstore.nnz = nnz;
	    Lstore.nsuper = col_to_sup[n];
	    Lstore.nzval = nzval;
	    Lstore.nzval_colbeg = nzval_colbeg;
	    Lstore.nzval_colend = nzval_colend;
	    Lstore.rowind = rowind;
	    Lstore.rowind_colbeg = rowind_colbeg;
	    Lstore.rowind_colend = rowind_colend;
	    Lstore.col_to_sup = col_to_sup;
	    Lstore.sup_to_colbeg = sup_to_colbeg;
	    Lstore.sup_to_colend = sup_to_colend;

	}


	/*
	 * Diagnostic print of column "jcol" in the U/L factor.
	 */
	public static
	void
	dprint_lu_col(int pnum, String msg, int pcol, int jcol, int w, int pivrow,
		      int xprune[], GlobalLU_t Glu)
	{
	    int     i, k, fsupc;
	    int     xsup[], supno[];
	    int     xlsub[], xlsub_end[], lsub[];
	    double  lusup[];
	    int     xlusup[], xlusup_end[];

	    xsup    = Glu.xsup;
	    supno   = Glu.supno;
	    lsub    = Glu.lsub;
	    xlsub   = Glu.xlsub;
	    xlsub_end = Glu.xlsub_end;
	    lusup   = Glu.lusup;
	    xlusup  = Glu.xlusup;
	    xlusup_end = Glu.xlusup_end;

	    printf("(%d)%s fstcol %d,col %d,w %d: pivrow %d, supno %d, xprune %d\n",
		   pnum, msg, pcol, jcol, w, pivrow, supno[jcol], xprune[jcol]);

	    printf("(%d)\tU-col: xusub %d - %d\n",
		   pnum, Glu.xusub[jcol], Glu.xusub_end[jcol]);
	    for (i = Glu.xusub[jcol]; i < Glu.xusub_end[jcol]; i++)
		printf("(%d)\t%d\t%8e\n", pnum, Glu.usub[i], Glu.ucol[i]);
	    fsupc = xsup[supno[jcol]];
	    k = xlusup[jcol];
	    printf("(%d)\tL-col in s-node: xlsub %d - %d, xlusup %d - %d\n",
		   pnum, xlsub[fsupc],xlsub_end[fsupc],xlusup[jcol],xlusup_end[jcol]);
	    for (i = xlsub[fsupc]; i < xlsub_end[fsupc]; ++i)
		printf("(%d)\t%d\t%.8e\n", pnum, lsub[i], lusup[k++]);

	    fflush(stdout);
	}

	/*Dan fix above printf 's*/

	/*
	 * Check whether vec[*] == 0. For the two vectors dense[*] and tempv[*],
	 * this invariant should be mantained before and after calling some
	 * numeric update routines, such as "panel_bmod" and "column_bmod".
	 */
	public static
	void
	dcheck_zero_vec(int pnum, String msg, int n, double vec[])
	{
	    int i, nonzero;

	    nonzero = FALSE;
	    for (i = 0; i < n; ++i) {
	        if (vec[i] != 0.0) {
	            printf("(%d) vec[%d] = %.10e; should be zero!\n",
	                   pnum, i, vec[i]);
	            nonzero = TRUE;
	        }
	    }
	    if ( nonzero != 0 ) {
		printf("(%d) %s\n", pnum, msg);
		SUPERLU_ABORT("Not a zero vector.");
	    }
	}


	public static
	void
	dGenXtrue(int n, int nrhs, double x[], int ldx)
	{
	    int  i, j;
	    for (j = 0; j < nrhs; ++j) {
		for (i = 0; i < n; ++i) {
	            x[i + j*ldx] = 1.0;/* + (double)(i+1.)/n;*/
	        }
	    }
	}

	/*
	 * Let rhs[i] = sum of i-th row of A, so the solution vector is all 1's
	 */
	public static
	void
	dFillRHS(trans_t trans, int nrhs, double x[], int ldx, SuperMatrix A, SuperMatrix B)
	{
	    DNformat Bstore;
	    double   rhs[];
	    double one = 1.0;
	    double zero = 0.0;
	    int      ldc;
	    char     trans_c;

	    Bstore = (DNformat) B.Store;
	    rhs    = Bstore.nzval;
	    ldc    = Bstore.lda;

	    if ( trans == NOTRANS ) trans_c = 'N';
	    else trans_c = 'T';

	    sp_dgemm(trans_c, A.nrow, nrhs, A.ncol, one, A,
		     x, ldx, zero, rhs, ldc);
	}

	/*
	 * Fills a double precision array with a given value.
	 */
	public static
	void
	dfill(double a[], int alen, double dval)
	{
	    int i;
	    for (i = 0; i < alen; i++) a[i] = dval;
	}



	/*
	 * Check the inf-norm of the error vector
	 */
	public static
	void dinf_norm_error(int nrhs, SuperMatrix X, double xtrue[])
	{
	    DNformat Xstore;
	    double err, xnorm;
	    double Xmat[], soln_work[];
	    int soln_work_offset;
	    int i, j;

	    Xstore = (DNformat) X.Store;
	    Xmat = Xstore.nzval;

	    for (j = 0; j < nrhs; j++) {
	      soln_work = Xmat;
	      soln_work_offset = j*Xstore.lda;
	      err = xnorm = 0.0;
	      for (i = 0; i < X.nrow; i++) {
	        err = SUPERLU_MAX(err, fabs(soln_work[soln_work_offset + i] - xtrue[i]));
	        xnorm = SUPERLU_MAX(xnorm, fabs(soln_work[soln_work_offset + i]));
	      }
	      err = err / xnorm;
	      printf("||X - Xtrue||/||X|| = %e\n", err);
	    }
	}



	/* Print performance of the code. */
	public static
	void
	dPrintPerf(SuperMatrix L, SuperMatrix U, superlu_memusage_t superlu_memusage,
		double rpg, double rcond, double ferr[],
		double berr[], String equed, Gstat_t Gstat)
	{
	    SCPformat Lstore;
	    NCPformat Ustore;
	    double   utime[];
	    float  ops[];

	    utime = Gstat.utime;
	    ops   = Gstat.ops;

	    if ( utime[FACT.ordinal()] != 0. )
		printf("Factor flops = %e\tMflops = %8.2f\n", ops[FACT.ordinal()],
		       ops[FACT.ordinal()]*1e-6/utime[FACT.ordinal()]);
	    printf("Identify relaxed snodes	= %8.2f\n", utime[RELAX.ordinal()]);
	    if ( utime[SOLVE.ordinal()] != 0. )
		printf("Solve flops = %.0f, Mflops = %8.2f\n", ops[SOLVE.ordinal()],
		       ops[SOLVE.ordinal()]*1e-6/utime[SOLVE.ordinal()]);

	    Lstore = (SCPformat) L.Store;
	    Ustore = (NCPformat) U.Store;
	    printf("\t#NZ in factor L = %d\n", Lstore.nnz);
	    printf("\t#NZ in factor U = %d\n", Ustore.nnz);
	    printf("\t#NZ in L+U = %d\n", Lstore.nnz + Ustore.nnz - L.ncol);

	    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
		   superlu_memusage.for_lu/1e6, superlu_memusage.total_needed/1e6,
		   superlu_memusage.expansions);

	    printf("\tFactor\tMflops\tSolve\tMflops\tEtree\tEquil\tRcond\tRefine\n");
	    printf("PERF:%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n",
		   utime[FACT.ordinal()], ops[FACT.ordinal()]*1e-6/utime[FACT.ordinal()],
		   utime[SOLVE.ordinal()], ops[SOLVE.ordinal()]*1e-6/utime[SOLVE.ordinal()],
		   utime[ETREE.ordinal()], utime[EQUIL.ordinal()], utime[RCOND.ordinal()],
		   utime[REFINE.ordinal()]);

	    printf("\tRpg\t\tRcond\t\tFerr\t\tBerr\t\tEquil?\n");
	    printf("NUM:\t%e\t%e\t%e\t%e\t%s\n",
		   rpg, rcond, ferr[0], berr[0], equed);

//	if (false) {
//
//	    printf("\tTRSV (total%%)\tGEMV (total%%)\tfloat_time%%\tmax_n\tmax_m\tmin_n\tmin_m\tavg_n\tavg_m\n");
//	    printf("BLAS:\t%.0f  %.2f\t%.0f  %.2f\t%.2f\t\t%d\t%d\t%d\t%d\t%.0f\t%.0f\n",
//		   ops[TRSV.ordinal()], ops[TRSV.ordinal()]/ops[FACT.ordinal()], ops[GEMV.ordinal()],
//		   ops[GEMV.ordinal()]/ops[FACT.ordinal()],
//		   utime[FLOAT.ordinal()]/utime[FACT.ordinal()],
//		   max_blas_n, max_gemv_m, min_blas_n, min_gemv_m,
//		   (float)sum_blas_n/num_blas, (float)sum_gemv_m/num_blas);
//	    printf("\tRCOND\tREFINE\tFERR\n");
//	    printf("SOLVES:\t%d\t%d\t%d\n", no_solves[RCOND.ordinal()],
//		   no_solves[REFINE.ordinal()], no_solves[FERR.ordinal()]);
//
//	    flops_dist_for_matlab();
//
//	}

	}


	public static
	int print_double_vec(String what, int n, int ind[], double vec[])
	{
		return print_double_vec(what, n, ind, 0, vec, 0);
	}

	public static
	int print_double_vec(String what, int n, int ind[], int ind_offset,
			double vec[], int vec_offset)
	{
	    int i;
	    printf("%s: n %d\n", what, n);
	    for (i = 0; i < n; ++i) printf("%d\t%f\n", ind[ind_offset+i], vec[vec_offset+i]);
	    return 0;
	}

}

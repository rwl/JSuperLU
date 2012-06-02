package gov.lbl.superlu.test;

import gov.lbl.superlu.Dlu_supermatrix.NCPformat;
import gov.lbl.superlu.Dlu_supermatrix.NCformat;
import gov.lbl.superlu.Dlu_supermatrix.SCPformat;
import gov.lbl.superlu.Dlu_supermatrix.SuperMatrix;

import gov.lbl.superlu.Dlu_slu_mt_util.superlu_memusage_t;
import gov.lbl.superlu.Dlu_slu_mt_util.trans_t;

import static gov.lbl.superlu.Dlu_sp_ienv.sp_ienv;

import static gov.lbl.superlu.Dlu_slu_mt_util.trans_t.NOTRANS;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_ABORT;

import static gov.lbl.superlu.Dlu_supermatrix.Stype_t.SLU_NC;
import static gov.lbl.superlu.Dlu_supermatrix.Stype_t.SLU_DN;
import static gov.lbl.superlu.Dlu_supermatrix.Dtype_t.SLU_D;
import static gov.lbl.superlu.Dlu_supermatrix.Mtype_t.SLU_GE;

import static gov.lbl.superlu.Dlu_pdutil.dCreate_CompCol_Matrix;
import static gov.lbl.superlu.Dlu_pdutil.dCreate_Dense_Matrix;
import static gov.lbl.superlu.Dlu_pdutil.dGenXtrue;
import static gov.lbl.superlu.Dlu_pdutil.dFillRHS;
import static gov.lbl.superlu.Dlu_pdutil.dinf_norm_error;

import static gov.lbl.superlu.Dlu_pdmemory.doubleMalloc;
import static gov.lbl.superlu.Dlu_pdmemory.superlu_dQuerySpace;

import static gov.lbl.superlu.Dlu_pmemory.intMalloc;

import static gov.lbl.superlu.Dlu.printf;

import static gov.lbl.superlu.Dlu_get_perm_c.get_perm_c;

import static gov.lbl.superlu.Dlu_pdgssv.pdgssv;


public class Dlu_pdlinsol {

	/* maximum number of processors to use. */
	public static final int nprocs = Runtime.getRuntime().availableProcessors();

	public static final int nrhs              = 1;
	public static final trans_t trans             = NOTRANS;
	public static final int n                 = 1000;
	public static final int b                 = 1;
	public static final int panel_size        = sp_ienv(1);
	public static final int relax             = sp_ienv(2);
	public static final int maxsup            = sp_ienv(3);

	static boolean HB = true;
	static boolean DEN = false;
	static boolean BAND = false;
	static boolean BD = false;

	public static void main(String[] args) {
		SuperMatrix   A = new SuperMatrix();
	    NCformat Astore;
	    double   a[];
	    int      asub[], xa[];
	    int      perm_r[]; /* row permutations from partial pivoting */
	    int      perm_c[]; /* column permutation vector */
	    SuperMatrix   L = new SuperMatrix();       /* factor L */
	    SCPformat Lstore;
	    SuperMatrix   U = new SuperMatrix();       /* factor U */
	    NCPformat Ustore;
	    SuperMatrix   B = new SuperMatrix();
	    int      ldx, info[], m, nnz;
	    info = new int[1];
	    int      permc_spec;
	    double   xact[], rhs[];
	    superlu_memusage_t   superlu_memusage;
//	    void   parse_command_line();


//	    parse_command_line(argc, argv, &nprocs, &n, &b, &panel_size,
//			       &relax, &maxsup);


//	if ( DEN ) {
//	    m = n;
//	    nnz = n * n;
//	    dband(n, n, nnz, &a, &asub, &xa);
//	} else if ( BAND ) {
//	    m = n;
//	    nnz = (2*b+1) * n;
//	    dband(n, b, nnz, &a, &asub, &xa);
//	} else if ( BD ) {
//	    nb = 5;
//	    bs = 200;
//	    m = n = bs * nb;
//	    nnz = bs * bs * nb;
//	    dblockdiag(nb, bs, nnz, &a, &asub, &xa);
//	} else if ( HB ) {
//	    dreadhb(&m, &n, &nnz, &a, &asub, &xa);
//	} else {
//	    dreadmt(&m, &n, &nnz, &a, &asub, &xa);
//	}

	    dCreate_CompCol_Matrix(A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
	    Astore = (NCformat) A.Store;
	    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore.nnz);

	    if ((rhs = doubleMalloc(m * nrhs)) == null) SUPERLU_ABORT("Malloc fails for rhs[].");
	    dCreate_Dense_Matrix(B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
	    xact = doubleMalloc(n * nrhs);
	    ldx = n;
	    dGenXtrue(n, nrhs, xact, ldx);
	    dFillRHS(trans, nrhs, xact, ldx, A, B);

	    if ((perm_r = intMalloc(m)) == null) SUPERLU_ABORT("Malloc fails for perm_r[].");
	    if ((perm_c = intMalloc(n)) == null) SUPERLU_ABORT("Malloc fails for perm_c[].");

	    /*
	     * Get column permutation vector perm_c[], according to permc_spec:
	     *   permc_spec = 0: natural ordering
	     *   permc_spec = 1: minimum degree ordering on structure of A'*A
	     *   permc_spec = 2: minimum degree ordering on structure of A'+A
	     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	     */
	    permc_spec = 1;
	    get_perm_c(permc_spec, A, perm_c);

	    pdgssv(nprocs, A, perm_c, perm_r, L, U, B, info);

	    if ( info[0] == 0 ) {
		dinf_norm_error(nrhs, B, xact); /* Inf. norm of the error */

		Lstore = (SCPformat) L.Store;
		Ustore = (NCPformat) U.Store;
	    	printf("#NZ in factor L = %d\n", Lstore.nnz);
	    	printf("#NZ in factor U = %d\n", Ustore.nnz);
	    	printf("#NZ in L+U = %d\n", Lstore.nnz + Ustore.nnz - L.ncol);

		superlu_dQuerySpace(nprocs, L, U, panel_size, superlu_memusage);
		printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
		       superlu_memusage.for_lu/1024/1024,
		       superlu_memusage.total_needed/1024/1024,
		       superlu_memusage.expansions);

	    }

	}

}

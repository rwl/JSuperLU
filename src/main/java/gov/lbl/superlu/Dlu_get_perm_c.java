package gov.lbl.superlu;

import org.netlib.util.intW;

import gov.lbl.superlu.Dlu_supermatrix.NCformat;
import gov.lbl.superlu.Dlu_supermatrix.SuperMatrix;

import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_MAX;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_ABORT;
import static gov.lbl.superlu.Dlu_slu_mt_util.FALSE;

import static gov.lbl.superlu.Dlu_superlu_timer.SuperLU_timer_;

import static gov.lbl.superlu.Dlu.printf;

import static edu.ufl.cise.colamd.tdouble.Dcolamd.COLAMD_STATS;
import static edu.ufl.cise.colamd.tdouble.Dcolamd.COLAMD_recommended;
import static edu.ufl.cise.colamd.tdouble.Dcolamd.COLAMD_set_defaults;
import static edu.ufl.cise.colamd.tdouble.Dcolamd.COLAMD_KNOBS;
import static edu.ufl.cise.colamd.tdouble.Dcolamd.colamd;

import static gov.lbl.superlu.mmd.Genmmd.genmmd;


public class Dlu_get_perm_c {

	public static
	void
	get_colamd(
		   final int m,  /* number of rows in matrix A. */
		   final int n,  /* number of columns in matrix A. */
		   final int nnz,/* number of nonzeros in matrix A. */
		   int colptr[],  /* column pointer of size n+1 for matrix A. */
		   int rowind[],  /* row indices of size nz for matrix A. */
		   int perm_c[]   /* out - the column permutation vector. */
		   )
	{
	    int Alen, A[], i, info, p[];
	    double knobs[];
		int[] stats = new int [COLAMD_STATS] ;

	    Alen = COLAMD_recommended(nnz, m, n);

	    if ( (knobs = new double[COLAMD_KNOBS]) == null )
	        SUPERLU_ABORT("Malloc fails for knobs");
	    COLAMD_set_defaults(knobs);

	    if (((A = (new int[Alen])) == null))
	        SUPERLU_ABORT("Malloc fails for A[]");
	    if ((p = (new int [n+1])) == null)
	        SUPERLU_ABORT("Malloc fails for p[]");
	    for (i = 0; i <= n; ++i) p[i] = colptr[i];
	    for (i = 0; i < nnz; ++i) A[i] = rowind[i];
	    info = colamd(m, n, Alen, A, p, knobs, stats);
	    if ( info == FALSE ) SUPERLU_ABORT("COLAMD failed");

	    for (i = 0; i < n; ++i) perm_c[p[i]] = i;
	}

	public static
	void
	getata(
	       final int m,      /* number of rows in matrix A. */
	       final int n,      /* number of columns in matrix A. */
	       final int nz,     /* number of nonzeros in matrix A */
	       int colptr[],      /* column pointer of size n+1 for matrix A. */
	       int rowind[],      /* row indices of size nz for matrix A. */
	       int atanz[],       /* out - on exit, returns the actual number of
	                            nonzeros in matrix A'*A. */
	       int ata_colptr[][], /* out - size n+1 */
	       int ata_rowind[][]  /* out - size *atanz */
	       )
	/*
	 * Purpose
	 * =======
	 *
	 * Form the structure of A'*A. A is an m-by-n matrix in column oriented
	 * format represented by (colptr, rowind). The output A'*A is in column
	 * oriented format (symmetrically, also row oriented), represented by
	 * (ata_colptr, ata_rowind).
	 *
	 * This routine is modified from GETATA routine by Tim Davis.
	 * The complexity of this algorithm is: SUM_{i=1,m} r(i)^2,
	 * i.e., the sum of the square of the row counts.
	 *
	 * Questions
	 * =========
	 *     o  Do I need to withhold the *dense* rows?
	 *     o  How do I know the number of nonzeros in A'*A?
	 *
	 */
	{
	    int i, j, k, col, num_nz, ti, trow;
	    int marker[], b_colptr[], b_rowind[];
	    int t_colptr[], t_rowind[]; /* a column oriented form of T = A' */

	    if ((marker = (new int[(SUPERLU_MAX(m,n)+1)])) == null)
		SUPERLU_ABORT("SUPERLU_MALLOC fails for marker[]");
	    if ( (t_colptr = (new int[m+1])) == null )
		SUPERLU_ABORT("SUPERLU_MALLOC t_colptr[]");
	    if ( (t_rowind = (new int[nz])) == null )
		SUPERLU_ABORT("SUPERLU_MALLOC fails for t_rowind[]");


	    /* Get counts of each column of T, and set up column pointers */
	    for (i = 0; i < m; ++i) marker[i] = 0;
	    for (j = 0; j < n; ++j) {
		for (i = colptr[j]; i < colptr[j+1]; ++i)
		    ++marker[rowind[i]];
	    }
	    t_colptr[0] = 0;
	    for (i = 0; i < m; ++i) {
		t_colptr[i+1] = t_colptr[i] + marker[i];
		marker[i] = t_colptr[i];
	    }

	    /* Transpose the matrix from A to T */
	    for (j = 0; j < n; ++j)
		for (i = colptr[j]; i < colptr[j+1]; ++i) {
		    col = rowind[i];
		    t_rowind[marker[col]] = j;
		    ++marker[col];
		}


	    /* ----------------------------------------------------------------
	       compute B = T * A, where column j of B is:

	       Struct (B_*j) =    UNION   ( Struct (T_*k) )
	                        A_kj != 0

	       do not include the diagonal entry

	       ( Partition A as: A = (A_*1, ..., A_*n)
	         Then B = T * A = (T * A_*1, ..., T * A_*n), where
	         T * A_*j = (T_*1, ..., T_*m) * A_*j.  )
	       ---------------------------------------------------------------- */

	    /* Zero the diagonal flag */
	    for (i = 0; i < n; ++i) marker[i] = -1;

	    /* First pass determines number of nonzeros in B */
	    num_nz = 0;
	    for (j = 0; j < n; ++j) {
		/* Flag the diagonal so it's not included in the B matrix */
		marker[j] = j;

		for (i = colptr[j]; i < colptr[j+1]; ++i) {
		    /* A_kj is nonzero, add pattern of column T_*k to B_*j */
		    k = rowind[i];
		    for (ti = t_colptr[k]; ti < t_colptr[k+1]; ++ti) {
			trow = t_rowind[ti];
			if ( marker[trow] != j ) {
			    marker[trow] = j;
			    num_nz++;
			}
		    }
		}
	    }
	    atanz[0] = num_nz;

	    /* Allocate storage for A'*A */
	    if ( (ata_colptr[0] = (new int[n+1])) == null )
		SUPERLU_ABORT("SUPERLU_MALLOC fails for ata_colptr[]");
	    if ( atanz[0] != 0 ) {
		if ( (ata_rowind[0] = (new int[atanz[0]])) == null )
		    SUPERLU_ABORT("SUPERLU_MALLOC fails for ata_rowind[]");
	    }
	    b_colptr = ata_colptr[0]; /* aliasing */
	    b_rowind = ata_rowind[0];

	    /* Zero the diagonal flag */
	    for (i = 0; i < n; ++i) marker[i] = -1;

	    /* Compute each column of B, one at a time */
	    num_nz = 0;
	    for (j = 0; j < n; ++j) {
		b_colptr[j] = num_nz;

		/* Flag the diagonal so it's not included in the B matrix */
		marker[j] = j;

		for (i = colptr[j]; i < colptr[j+1]; ++i) {
		    /* A_kj is nonzero, add pattern of column T_*k to B_*j */
		    k = rowind[i];
		    for (ti = t_colptr[k]; ti < t_colptr[k+1]; ++ti) {
			trow = t_rowind[ti];
			if ( marker[trow] != j ) {
			    marker[trow] = j;
			    b_rowind[num_nz++] = trow;
			}
		    }
		}
	    }
	    b_colptr[n] = num_nz;
	}


	public static
	void
	at_plus_a(
		  final int n,      /* number of columns in matrix A. */
		  final int nz,     /* number of nonzeros in matrix A */
		  int colptr[],      /* column pointer of size n+1 for matrix A. */
		  int rowind[],      /* row indices of size nz for matrix A. */
		  int bnz[],         /* out - on exit, returns the actual number of
	                               nonzeros in matrix A'*A. */
		  int b_colptr[][],   /* out - size n+1 */
		  int b_rowind[][]    /* out - size *bnz */
		  )
	{
	/*
	 * Purpose
	 * =======
	 *
	 * Form the structure of A'+A. A is an n-by-n matrix in column oriented
	 * format represented by (colptr, rowind). The output A'+A is in column
	 * oriented format (symmetrically, also row oriented), represented by
	 * (b_colptr, b_rowind).
	 *
	 */
	    int i, j, k, col, num_nz;
	    int t_colptr[], t_rowind[]; /* a column oriented form of T = A' */
	    int marker[];

	    if ( (marker = (new int[n])) == null )
		SUPERLU_ABORT("SUPERLU_MALLOC fails for marker[]");
	    if ( (t_colptr = (new int[n+1])) == null )
		SUPERLU_ABORT("SUPERLU_MALLOC fails for t_colptr[]");
	    if ( (t_rowind = (new int[nz])) == null )
		SUPERLU_ABORT("SUPERLU_MALLOC fails t_rowind[]");


	    /* Get counts of each column of T, and set up column pointers */
	    for (i = 0; i < n; ++i) marker[i] = 0;
	    for (j = 0; j < n; ++j) {
		for (i = colptr[j]; i < colptr[j+1]; ++i)
		    ++marker[rowind[i]];
	    }
	    t_colptr[0] = 0;
	    for (i = 0; i < n; ++i) {
		t_colptr[i+1] = t_colptr[i] + marker[i];
		marker[i] = t_colptr[i];
	    }

	    /* Transpose the matrix from A to T */
	    for (j = 0; j < n; ++j)
		for (i = colptr[j]; i < colptr[j+1]; ++i) {
		    col = rowind[i];
		    t_rowind[marker[col]] = j;
		    ++marker[col];
		}


	    /* ----------------------------------------------------------------
	       compute B = A + T, where column j of B is:

	       Struct (B_*j) = Struct (A_*k) UNION Struct (T_*k)

	       do not include the diagonal entry
	       ---------------------------------------------------------------- */

	    /* Zero the diagonal flag */
	    for (i = 0; i < n; ++i) marker[i] = -1;

	    /* First pass determines number of nonzeros in B */
	    num_nz = 0;
	    for (j = 0; j < n; ++j) {
		/* Flag the diagonal so it's not included in the B matrix */
		marker[j] = j;

		/* Add pattern of column A_*k to B_*j */
		for (i = colptr[j]; i < colptr[j+1]; ++i) {
		    k = rowind[i];
		    if ( marker[k] != j ) {
			marker[k] = j;
			++num_nz;
		    }
		}

		/* Add pattern of column T_*k to B_*j */
		for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
		    k = t_rowind[i];
		    if ( marker[k] != j ) {
			marker[k] = j;
			++num_nz;
		    }
		}
	    }
	    bnz[0] = num_nz;

	    /* Allocate storage for A+A' */
	    if ( (b_colptr[0] = (new int[n+1])) == null )
		SUPERLU_ABORT("SUPERLU_MALLOC fails for b_colptr[]");
	    if ( bnz[0] != 0 ) {
	      if ( (b_rowind[0] = (new int[bnz[0]])) == null )
		SUPERLU_ABORT("SUPERLU_MALLOC fails for b_rowind[]");
	    }

	    /* Zero the diagonal flag */
	    for (i = 0; i < n; ++i) marker[i] = -1;

	    /* Compute each column of B, one at a time */
	    num_nz = 0;
	    for (j = 0; j < n; ++j) {
		(b_colptr[0])[j] = num_nz;

		/* Flag the diagonal so it's not included in the B matrix */
		marker[j] = j;

		/* Add pattern of column A_*k to B_*j */
		for (i = colptr[j]; i < colptr[j+1]; ++i) {
		    k = rowind[i];
		    if ( marker[k] != j ) {
			marker[k] = j;
			(b_rowind[0])[num_nz++] = k;
		    }
		}

		/* Add pattern of column T_*k to B_*j */
		for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
		    k = t_rowind[i];
		    if ( marker[k] != j ) {
			marker[k] = j;
			(b_rowind[0])[num_nz++] = k;
		    }
		}
	    }
	    (b_colptr[0])[n] = num_nz;
	}

	public static
	void
	get_perm_c(int ispec, SuperMatrix A, int perm_c[])
	/*
	 * Purpose
	 * =======
	 *
	 * GET_PERM_C obtains a permutation matrix Pc, by applying the multiple
	 * minimum degree ordering code by Joseph Liu to matrix A'*A or A+A'.
	 * or using approximate minimum degree column ordering by Davis et. al.
	 * The LU factorization of A*Pc tends to have less fill than the LU
	 * factorization of A.
	 *
	 * Arguments
	 * =========
	 *
	 * ispec   (input) int
	 *         Specifies the type of column ordering to reduce fill:
	 *         = 1: minimum degree on the structure of A^T * A
	 *         = 2: minimum degree on the structure of A^T + A
	 *         = 3: approximate minimum degree for unsymmetric matrices
	 *         If ispec == 0, the natural ordering (i.e., Pc = I) is returned.
	 *
	 * A       (input) SuperMatrix*
	 *         Matrix A in A*X=B, of dimension (A.nrow, A.ncol). The number
	 *         of the linear equations is A.nrow. Currently, the type of A
	 *         can be: Stype = NC; Dtype = _D; Mtype = GE. In the future,
	 *         more general A can be handled.
	 *
	 * perm_c  (output) int*
	 *	   Column permutation vector of size A.ncol, which defines the
	 *         permutation matrix Pc; perm_c[i] = j means column i of A is
	 *         in position j in A*Pc.
	 *
	 */
	{
	    NCformat Astore = (NCformat) A.Store;
	    int m, n, bnz[], b_colptr[][], i;
	    bnz = new int[1];
	    b_colptr = new int[1][];
	    int delta, maxint, invp[];
	    intW nofsub = new intW(0);
	    int b_rowind[][], dhead[], qsize[], llist[], marker[];
	    b_rowind = new int[1][];
	    double t;

	    m = A.nrow;
	    n = A.ncol;

	    t = SuperLU_timer_();
	    switch ( ispec ) {
	        case 0: /* Natural ordering */
		      for (i = 0; i < n; ++i) perm_c[i] = i;
		      printf("Use natural column ordering.\n");
		      return;
	        case 1: /* Minimum degree ordering on A'*A */
		      getata(m, n, Astore.nnz, Astore.colptr, Astore.rowind,
			     bnz, b_colptr, b_rowind);
		      printf("Use minimum degree ordering on A'*A.\n");
		      t = SuperLU_timer_() - t;
		      /*printf("Form A'*A time = %8.3f\n", t);*/
		      break;
	        case 2: /* Minimum degree ordering on A'+A */
		      if ( m != n ) SUPERLU_ABORT("Matrix is not square");
		      at_plus_a(n, Astore.nnz, Astore.colptr, Astore.rowind,
				bnz, b_colptr, b_rowind);
		      printf("Use minimum degree ordering on A'+A.\n");
		      t = SuperLU_timer_() - t;
		      /*printf("Form A'+A time = %8.3f\n", t);*/
		      break;
	        case 3: /* Approximate minimum degree column ordering. */
		      get_colamd(m, n, Astore.nnz, Astore.colptr, Astore.rowind,
				 perm_c);
		      printf(".. Use approximate minimum degree column ordering.\n");
		      return;
	        default:
		      SUPERLU_ABORT("Invalid ISPEC");
	    }

	    if ( bnz[0] != 0 ) {
		t = SuperLU_timer_();

		/* Initialize and allocate storage for GENMMD. */
		delta = 0; /* DELTA is a parameter to allow the choice of nodes
			      whose degree <= min-degree + DELTA. */
		maxint = 2147483647; /* 2**31 - 1 */
		invp = (new int[n+delta]);
		if ( invp == null ) SUPERLU_ABORT("SUPERLU_MALLOC fails for invp.");
		dhead = (new int [n+delta]);
		if ( dhead == null ) SUPERLU_ABORT("SUPERLU_MALLOC fails for dhead.");
		qsize = (new int [n+delta]);
		if ( qsize == null ) SUPERLU_ABORT("SUPERLU_MALLOC fails for qsize.");
		llist = (new int [n]);
		if ( llist == null ) SUPERLU_ABORT("SUPERLU_MALLOC fails for llist.");
		marker = (new int [n]);
		if ( marker == null ) SUPERLU_ABORT("SUPERLU_MALLOC fails for marker.");

		/* Transform adjacency list into 1-based indexing required by GENMMD.*/
		for (i = 0; i <= n; ++i) ++b_colptr[0][i];
		for (i = 0; i < bnz[0]; ++i) ++b_rowind[0][i];

		genmmd(n, b_colptr[0], 0, b_rowind[0], 0, perm_c, 0, invp, 0, delta, dhead, 0,
			qsize, 0, llist, 0, marker, 0, maxint, nofsub);

		/* Transform perm_c into 0-based indexing. */
		for (i = 0; i < n; ++i) --perm_c[i];

		t = SuperLU_timer_() - t;
		/*  printf("call GENMMD time = %8.3f\n", t);*/

	    } else { /* Empty adjacency structure */
		for (i = 0; i < n; ++i) perm_c[i] = i;
	    }

	}
}

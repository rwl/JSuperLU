/*! @file Dlu_supermatrix.java
 * \brief Defines matrix types
 */
package lu.jsuper;


public class Dlu_supermatrix {

	/********************************************
	 * The matrix types are defined as follows. *
	 ********************************************/
	public enum Stype_t {
	    SLU_NC,    /* column-wise, no supernode */
	    SLU_NCP,   /* column-wise, column-permuted, no supernode
	                  (The consecutive columns of nonzeros, after permutation,
			   may not be stored  contiguously.) */
	    SLU_NR,    /* row-wize, no supernode */
	    SLU_SC,    /* column-wise, supernode */
	    SLU_SCP,   /* supernode, column-wise, permuted */
	    SLU_SR,    /* row-wise, supernode */
	    SLU_DN,     /* Fortran style column-wise storage for dense matrix */
	    SLU_NR_loc; /* distributed compressed row format  */
	}

	public enum Dtype_t {
	    SLU_S,     /* single */
	    SLU_D,     /* double */
	    SLU_C,     /* single complex */
	    SLU_Z;     /* double complex */
	}

	public enum Mtype_t {
	    SLU_GE,    /* general */
	    SLU_TRLU,  /* lower triangular, unit diagonal */
	    SLU_TRUU,  /* upper triangular, unit diagonal */
	    SLU_TRL,   /* lower triangular */
	    SLU_TRU,   /* upper triangular */
	    SLU_SYL,   /* symmetric, store lower half */
	    SLU_SYU,   /* symmetric, store upper half */
	    SLU_HEL,   /* Hermitian, store lower half */
	    SLU_HEU;   /* Hermitian, store upper half */
	}

	public static class SuperMatrix {
		public 	Stype_t Stype; /* Storage type: interprets the storage structure
			pointed to by *Store. */
		public Dtype_t Dtype; /* Data type. */
		public Mtype_t Mtype; /* Matrix type: describes the mathematical property of
			the matrix. */
		public int nrow;   /* number of rows */
		public int ncol;   /* number of columns */
		Object Store;   /* pointer to the actual storage of the matrix */
	}

	/***********************************************
	 * The storage schemes are defined as follows. *
	 ***********************************************/

	/* Stype == SLU_NC (Also known as Harwell-Boeing sparse matrix format) */
	public static class NCformat {
	    int    nnz;	     /* number of nonzeros in the matrix */
	    double nzval[];  /* pointer to array of nonzero values, packed by column */
	    int    rowind[]; /* pointer to array of row indices of the nonzeros */
	    int    colptr[]; /* pointer to array of beginning of columns in nzval[]
			       and rowind[]  */
	            /* Note:
			       Zero-based indexing is used;
			       colptr[] has ncol+1 entries, the last one pointing
			       beyond the last column, so that colptr[ncol] = nnz. */
	}

	/* Stype == SLU_NR */
	public static class NRformat {
	    int    nnz;	    /* number of nonzeros in the matrix */
	    double nzval[]; /* pointer to array of nonzero values, packed by raw */
	    int  colind[];  /* pointer to array of columns indices of the nonzeros */
	    int  rowptr[];  /* pointer to array of beginning of rows in nzval[]
			       and colind[]  */
	                    /* Note:
			       Zero-based indexing is used;
			       rowptr[] has nrow+1 entries, the last one pointing
			       beyond the last row, so that rowptr[nrow] = nnz. */
	}

	/* Stype == SLU_DN */
	public static class DNformat {
	    int lda;    /* leading dimension */
	    double nzval[];  /* array of size lda*ncol to represent a dense matrix */
	}

}

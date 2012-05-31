package gov.lbl.superlu;

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
	    SLU_NR_loc  /* distributed compressed row format  */
	}

	public enum Dtype_t {
	    SLU_S,     /* single */
	    SLU_D,     /* double */
	    SLU_C,     /* single complex */
	    SLU_Z      /* double complex */
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
	    SLU_HEU    /* Hermitian, store upper half */
	}

	public static class SuperMatrix {
		Stype_t Stype; /* Storage type: interprets the storage structure
			   	  pointed to by *Store. */
		Dtype_t Dtype; /* Data type. */
		Mtype_t Mtype; /* Matrix type: describes the mathematical property of
				  the matrix. */
		int    nrow;   /* number of rows */
		int    ncol;   /* number of columns */
		Object Store;   /* pointer to the actual storage of the matrix */
	}

	/***********************************************
	 * The storage schemes are defined as follows. *
	 ***********************************************/

	/* Stype == SLU_NC (Also known as Harwell-Boeing sparse matrix format) */
	public static class NCformat{
	    int    nnz;	    /* number of nonzeros in the matrix */
	    double nzval[]; /* pointer to array of nonzero values, packed by column */
	    int    rowind[];/* pointer to array of row indices of the nonzeros */
	    int    colptr[];/* pointer to array of beginning of columns in nzval[]
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
	    int    colind[];/* pointer to array of columns indices of the nonzeros */
	    int    rowptr[];/* pointer to array of beginning of rows in nzval[]
			       and colind[]  */
	                    /* Note:
			       Zero-based indexing is used;
			       rowptr[] has nrow+1 entries, the last one pointing
			       beyond the last row, so that rowptr[nrow] = nnz. */
	}

	/* Stype == SLU_SC */
	public static class SCformat {
	  int    nnz;	     /* number of nonzeros in the matrix */
	  int    nsuper;     /* number of supernodes, minus 1 */
	  double nzval[];    /* pointer to array of nonzero values, packed by column */
	  int nzval_colptr[];/* pointer to array of beginning of columns in nzval[] */
	  int rowind[];      /* pointer to array of compressed row indices of
				rectangular supernodes */
	  int rowind_colptr[];/* pointer to array of beginning of columns in rowind[] */
	  int col_to_sup[];   /* col_to_sup[j] is the supernode number to which column
				j belongs; mapping from column to supernode number. */
	  int sup_to_col[];   /* sup_to_col[s] points to the start of the s-th
				supernode; mapping from supernode number to column.
			        e.g.: col_to_sup: 0 1 2 2 3 3 3 4 4 4 4 4 4 (ncol=12)
			              sup_to_col: 0 1 2 4 7 12           (nsuper=4) */
	                     /* Note:
			        Zero-based indexing is used;
			        nzval_colptr[], rowind_colptr[], col_to_sup and
			        sup_to_col[] have ncol+1 entries, the last one
			        pointing beyond the last column.
			        For col_to_sup[], only the first ncol entries are
			        defined. For sup_to_col[], only the first nsuper+2
			        entries are defined. */
	}

	/* Stype == SLU_SCP */
	public static class SCPformat {
	  int    nnz;	     /* number of nonzeros in the matrix */
	  int    nsuper;     /* number of supernodes */
	  double nzval[];    /* pointer to array of nonzero values, packed by column */
	  int    nzval_colbeg[];/* nzval_colbeg[j] points to beginning of column j
				  in nzval[] */
	  int    nzval_colend[];/* nzval_colend[j] points to one past the last element
				  of column j in nzval[] */
	  int    rowind[];      /* pointer to array of compressed row indices of
				  rectangular supernodes */
	  int   rowind_colbeg[];/* rowind_colbeg[j] points to beginning of column j
				  in rowind[] */
	  int   rowind_colend[];/* rowind_colend[j] points to one past the last element
				  of column j in rowind[] */
	  int   col_to_sup[];   /* col_to_sup[j] is the supernode number to which column
				  j belongs; mapping from column to supernode. */
	  int   sup_to_colbeg[]; /* sup_to_colbeg[s] points to the start of the s-th
				   supernode; mapping from supernode to column.*/
	  int   sup_to_colend[]; /* sup_to_colend[s] points to one past the end of the
				   s-th supernode; mapping from supernode number to
				   column.
			        e.g.: col_to_sup: 0 1 2 2 3 3 3 4 4 4 4 4 4 (ncol=12)
			              sup_to_colbeg: 0 1 2 4 7              (nsuper=4)
				      sup_to_colend: 1 2 4 7 12                    */
	                     /* Note:
			        Zero-based indexing is used;
			        nzval_colptr[], rowind_colptr[], col_to_sup and
			        sup_to_col[] have ncol+1 entries, the last one
			        pointing beyond the last column.         */
	}

	/* Stype == SLU_NCP */
	public static class NCPformat {
	    int   nnz;	  /* number of nonzeros in the matrix */
	    double nzval[];  /* pointer to array of nonzero values, packed by column */
	    int   rowind[];/* pointer to array of row indices of the nonzeros */
			  /* Note: nzval[]/rowind[] always have the same length */
	    int colbeg[];/* colbeg[j] points to the beginning of column j in nzval[]
	                     and rowind[]  */
	    int colend[];/* colend[j] points to one past the last element of column
			     j in nzval[] and rowind[]  */
			  /* Note:
			     Zero-based indexing is used;
			     The consecutive columns of the nonzeros may not be
			     contiguous in storage, because the matrix has been
			     postmultiplied by a column permutation matrix. */
	}

	/* Stype == SLU_DN */
	public static class DNformat {
	    int   lda;    /* leading dimension */
	    double nzval[];  /* array of size lda*ncol to represent a dense matrix */
	}

	/* Stype == SLU_NR_loc (Distributed Compressed Row Format) */
	public static class NRformat_loc {
	    int   nnz_loc;   /* number of nonzeros in the local submatrix */
	    int   m_loc;     /* number of rows local to this processor */
	    int   fst_row;   /* global index of the first row */
	    double nzval[];  /* pointer to array of nonzero values, packed by row */
	    int rowptr[];    /* pointer to array of beginning of rows in nzval[]
				and colind[]  */
	    int   colind[];   /* pointer to array of column indices of the nonzeros */
	                     /* Note:
				Zero-based indexing is used;
				rowptr[] has n_loc + 1 entries, the last one pointing
				beyond the last row, so that rowptr[n_loc] = nnz_loc.*/
	}

}

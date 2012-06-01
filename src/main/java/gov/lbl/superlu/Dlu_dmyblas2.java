/*! @file Dlu_dmyblas2.java
 * \brief Level 2 Blas operations
 *
 * <pre>
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 * </pre>
 * Purpose:
 *     Level 2 BLAS operations: solves and matvec, written in C.
 * Note:
 *     This is only used when the system lacks an efficient BLAS library.
 * </pre>
 */
package gov.lbl.superlu;

public class Dlu_dmyblas2 {

	/**! \brief Solves a dense UNIT lower triangular system
	 *
	 *  The unit lower
	 * triangular matrix is stored in a 2D array M(1:nrow,1:ncol).
	 * The solution will be returned in the rhs vector.
	 */
	public static void dlsolve(int ldm, int ncol, double M[], int M_offset,
			double rhs[], int rhs_offset) {
	    int k;
	    double x0, x1, x2, x3, x4, x5, x6, x7;
	    double[] M0;
	    int M0_offset;
	    double[] Mki0, Mki1, Mki2, Mki3, Mki4, Mki5, Mki6, Mki7;
	    int Mki0_offset, Mki1_offset, Mki2_offset, Mki3_offset, Mki4_offset,
	    Mki5_offset, Mki6_offset, Mki7_offset;
	    int firstcol = 0;

	    M0 = M;
	    M0_offset = M_offset;

	    while ( firstcol < ncol - 7 ) { /* Do 8 columns */
		  Mki0 = M0;
		  Mki0_offset = M0_offset + 1;
		  Mki1 = Mki0;
		  Mki1_offset = Mki0_offset + ldm + 1;
		  Mki2 = Mki1;
		  Mki2_offset = Mki1_offset + ldm + 1;
		  Mki3 = Mki2;
		  Mki3_offset = Mki2_offset + ldm + 1;
		  Mki4 = Mki3;
		  Mki4_offset = Mki3_offset + ldm + 1;
		  Mki5 = Mki4;
		  Mki5_offset = Mki4_offset + ldm + 1;
		  Mki6 = Mki5;
		  Mki6_offset = Mki5_offset + ldm + 1;
		  Mki7 = Mki6;
		  Mki7_offset = Mki6_offset + ldm + 1;

		  x0 = rhs[rhs_offset+firstcol];
		  x1 = rhs[rhs_offset+firstcol+1] - x0 * Mki0[Mki0_offset];
		  Mki0_offset++;
		  x2 = rhs[rhs_offset+firstcol+2] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset];
		  Mki0_offset++;
		  Mki1_offset++;
		  x3 = rhs[rhs_offset+firstcol+3] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
				  - x2 * Mki2[Mki2_offset];
		  Mki0_offset++;
		  Mki1_offset++;
		  Mki2_offset++;
		  x4 = rhs[rhs_offset+firstcol+4] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
				  - x2 * Mki2[Mki2_offset] - x3 * Mki3[Mki3_offset];
		  Mki0_offset++;
		  Mki1_offset++;
		  Mki2_offset++;
		  Mki3_offset++;
		  x5 = rhs[rhs_offset+firstcol+5] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
				  - x2 * Mki2[Mki2_offset] - x3 * Mki3[Mki3_offset]
				  - x4 * Mki4[Mki4_offset];
		  Mki0_offset++;
		  Mki1_offset++;
		  Mki2_offset++;
		  Mki3_offset++;
		  Mki4_offset++;
		  x6 = rhs[rhs_offset+firstcol+6] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
				  - x2 * Mki2[Mki2_offset] - x3 * Mki3[Mki3_offset]
				  - x4 * Mki4[Mki4_offset] - x5 * Mki5[Mki5_offset];
		  Mki0_offset++;
		  Mki1_offset++;
		  Mki2_offset++;
		  Mki3_offset++;
		  Mki4_offset++;
		  Mki5_offset++;
		  x7 = rhs[rhs_offset+firstcol+7] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
				  - x2 * Mki2[Mki2_offset] - x3 * Mki3[Mki3_offset]
				  - x4 * Mki4[Mki4_offset] - x5 * Mki5[Mki5_offset]
				  - x6 * Mki6[Mki6_offset];
		  Mki0_offset++;
		  Mki1_offset++;
		  Mki2_offset++;
		  Mki3_offset++;
		  Mki4_offset++;
		  Mki5_offset++;
		  Mki6_offset++;

		  rhs[rhs_offset+ ++firstcol] = x1;
		  rhs[rhs_offset+ ++firstcol] = x2;
		  rhs[rhs_offset+ ++firstcol] = x3;
		  rhs[rhs_offset+ ++firstcol] = x4;
		  rhs[rhs_offset+ ++firstcol] = x5;
		  rhs[rhs_offset+ ++firstcol] = x6;
		  rhs[rhs_offset+ ++firstcol] = x7;
		  ++firstcol;

	      for (k = firstcol; k < ncol; k++) {
	    	  rhs[rhs_offset+k] = rhs[rhs_offset+k] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
		                - x2 * Mki2[Mki2_offset] - x3 * Mki3[Mki3_offset]
	                    - x4 * Mki4[Mki4_offset] - x5 * Mki5[Mki5_offset]
	                    - x6 * Mki6[Mki6_offset] - x7 * Mki7[Mki7_offset];
	    	  Mki0_offset++;
	    	  Mki1_offset++;
	    	  Mki2_offset++;
	    	  Mki3_offset++;
	    	  Mki4_offset++;
	    	  Mki5_offset++;
	    	  Mki6_offset++;
	    	  Mki7_offset++;
	      }

	      M0_offset += 8 * ldm + 8;
	    }

	    while ( firstcol < ncol - 3 ) { /* Do 4 columns */
		  Mki0 = M0;
		  Mki0_offset = M0_offset + 1;
		  Mki1 = Mki0;
		  Mki1_offset = Mki0_offset + ldm + 1;
		  Mki2 = Mki1;
		  Mki2_offset = Mki1_offset + ldm + 1;
		  Mki3 = Mki2;
		  Mki3_offset = Mki2_offset + ldm + 1;

	      x0 = rhs[rhs_offset+firstcol];
	      x1 = rhs[rhs_offset+firstcol+1] - x0 * Mki0[Mki0_offset];
	      Mki0_offset++;
	      x2 = rhs[rhs_offset+firstcol+2] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset];
	      Mki0_offset++;
	      Mki1_offset++;
	      x3 = rhs[rhs_offset+firstcol+3] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
	    		  - x2 * Mki2[Mki2_offset];
	      Mki0_offset++;
	      Mki1_offset++;
	      Mki2_offset++;

	      rhs[rhs_offset+ ++firstcol] = x1;
	      rhs[rhs_offset+ ++firstcol] = x2;
	      rhs[rhs_offset+ ++firstcol] = x3;
	      ++firstcol;

	      for (k = firstcol; k < ncol; k++) {
	    	  rhs[rhs_offset+k] = rhs[rhs_offset+k] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
		                - x2 * Mki2[Mki2_offset] - x3 * Mki3[Mki3_offset];
	    	  Mki0_offset++;
	    	  Mki1_offset++;
	    	  Mki2_offset++;
	    	  Mki3_offset++;
	      }

	      M0_offset += 4 * ldm + 4;
	    }

	    if ( firstcol < ncol - 1 ) { /* Do 2 columns */
	      Mki0 = M0;
		  Mki0_offset = M0_offset + 1;
		  Mki1 = Mki0;
		  Mki1_offset = Mki0_offset + ldm + 1;

	      x0 = rhs[rhs_offset+firstcol];
	      x1 = rhs[rhs_offset+firstcol+1] - x0 * Mki0[Mki0_offset];
	      Mki0_offset++;

	      rhs[rhs_offset+ ++firstcol] = x1;
	      ++firstcol;

	      for (k = firstcol; k < ncol; k++) {
	    	  rhs[rhs_offset+k] = rhs[rhs_offset+k] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset];
	    	  Mki0_offset++;
	    	  Mki1_offset++;
	      }

	    }
	}

	/*! \brief Solves a dense upper triangular system
	 *
	 * The upper triangular matrix is
	 * stored in a 2-dim array M(1:ldm,1:ncol). The solution will be returned
	 * in the rhs vector.
	 */
	public static void dusolve(int ldm, int ncol, double M[], int M_offset,
			double rhs[], int rhs_offset)
	{
	    double xj;
	    int jcol, j, irow;

	    jcol = ncol - 1;

	    for (j = 0; j < ncol; j++) {

		xj = rhs[rhs_offset + jcol] / M[M_offset + jcol + jcol*ldm]; 		/* M(jcol, jcol) */
		rhs[rhs_offset + jcol] = xj;

		for (irow = 0; irow < jcol; irow++)
		    rhs[rhs_offset + irow] -= xj * M[M_offset + irow + jcol*ldm];	/* M(irow, jcol) */

		jcol--;

	    }
	}

	/*! \brief Performs a dense matrix-vector multiply: Mxvec = Mxvec + M * vec.
	 *
	 * The input matrix is M(1:nrow,1:ncol); The product is returned in Mxvec[].
	 */
	public static void dmatvec(
		int ldm,	/* in -- leading dimension of M */
		int nrow,	/* in */
		int ncol,	/* in */
		double M[],	/* in */
		int M_offset,
		double vec[],	/* in */
		int vec_offset,
		double Mxvec[],	/* in/out */
		int Mxvec_offset
		)
	{
	    double vi0, vi1, vi2, vi3, vi4, vi5, vi6, vi7;
	    double M0[];
	    int M0_offset;
	    double[] Mki0, Mki1, Mki2, Mki3, Mki4, Mki5, Mki6, Mki7;
	    int Mki0_offset, Mki1_offset, Mki2_offset, Mki3_offset, Mki4_offset,
	    Mki5_offset, Mki6_offset, Mki7_offset;
	    int firstcol = 0;
	    int k;

	    M0 = M;
	    M0_offset = M_offset;

	    while ( firstcol < ncol - 7 ) {	/* Do 8 columns */

		Mki0 = M0;
		Mki0_offset = M0_offset;
		Mki1 = Mki0;
		Mki1_offset = Mki0_offset + ldm;
	    Mki2 = Mki1;
	    Mki2_offset = Mki1_offset + ldm;
	    Mki3 = Mki2;
	    Mki3_offset = Mki2_offset + ldm;
		Mki4 = Mki3;
		Mki4_offset = Mki3_offset + ldm;
		Mki5 = Mki4;
		Mki5_offset = Mki4_offset + ldm;
		Mki6 = Mki5;
		Mki6_offset = Mki5_offset + ldm;
		Mki7 = Mki6;
		Mki7_offset = Mki6_offset + ldm;

		vi0 = vec[vec_offset+firstcol++];
		vi1 = vec[vec_offset+firstcol++];
		vi2 = vec[vec_offset+firstcol++];
		vi3 = vec[vec_offset+firstcol++];
		vi4 = vec[vec_offset+firstcol++];
		vi5 = vec[vec_offset+firstcol++];
		vi6 = vec[vec_offset+firstcol++];
		vi7 = vec[vec_offset+firstcol++];

		for (k = 0; k < nrow; k++) {
		    Mxvec[Mxvec_offset+k] += vi0 * Mki0[Mki0_offset] + vi1 * Mki1[Mki1_offset]
			      + vi2 * Mki2[Mki2_offset] + vi3 * Mki3[Mki3_offset]
			      + vi4 * Mki4[Mki4_offset] + vi5 * Mki5[Mki5_offset]
			      + vi6 * Mki6[Mki6_offset] + vi7 * Mki7[Mki7_offset];
		    Mki0_offset++;
		    Mki1_offset++;
		    Mki2_offset++;
		    Mki3_offset++;
		    Mki4_offset++;
		    Mki5_offset++;
		    Mki6_offset++;
		    Mki7_offset++;
		}

		M0_offset += 8 * ldm;
	    }

	    while ( firstcol < ncol - 3 ) {	/* Do 4 columns */

		Mki0 = M0;
		Mki0_offset = M0_offset;
		Mki1 = Mki0;
		Mki1_offset = Mki0_offset + ldm;
		Mki2 = Mki1;
		Mki2_offset = Mki1_offset + ldm;
		Mki3 = Mki2;
		Mki3_offset = Mki2_offset + ldm;

		vi0 = vec[vec_offset+firstcol++];
		vi1 = vec[vec_offset+firstcol++];
		vi2 = vec[vec_offset+firstcol++];
		vi3 = vec[vec_offset+firstcol++];
		for (k = 0; k < nrow; k++) {
		    Mxvec[Mxvec_offset+k] += vi0 * Mki0[Mki0_offset] + vi1 * Mki1[Mki1_offset]
			      + vi2 * Mki2[Mki2_offset] + vi3 * Mki3[Mki3_offset] ;
		    Mki0_offset++;
		    Mki1_offset++;
		    Mki2_offset++;
		    Mki3_offset++;
		}

		M0_offset += 4 * ldm;
	    }

	    while ( firstcol < ncol ) {		/* Do 1 column */

	 	Mki0 = M0;
	 	Mki0_offset = M0_offset;
		vi0 = vec[vec_offset+firstcol++];
		for (k = 0; k < nrow; k++) {
		    Mxvec[Mxvec_offset+k] += vi0 * Mki0[Mki0_offset];
		    Mki0_offset++;
		}

		M0_offset += ldm;
	    }

	}

	/*
	 * Performs dense matrix-vector multiply with 2 vectors:
	 *        y0 = y0 + A * x0
	 *        y1 = y1 + A * x1
	 */
	static
	void dmatvec2 (
	               int lda,     /* leading dimension of A */
	               int m,
	               int n,
	               double A[],   /* in - size m-by-n */
	               int A_offset,
	               double x0[],  /* in - size n-by-1 */
	               int x0_offset,
	               double x1[],  /* in - size n-by-1 */
	               int x1_offset,
	               double y0[],  /* out - size n-by-1 */
	               int y0_offset,
	               double y1[],   /* out - size n-by-1 */
	               int y1_offset
	               )

	{
	    double v00, v10, v20, v30, v40, v50, v60, v70,
	                    v01, v11, v21, v31, v41, v51, v61, v71;
	    double t0, t1, t2, t3, t4, t5, t6, t7;
	    double f0, f1;
	    double[] Mki0, Mki1, Mki2, Mki3, Mki4, Mki5, Mki6, Mki7;
	    int Mki0_offset, Mki1_offset, Mki2_offset, Mki3_offset,
	    Mki4_offset, Mki5_offset, Mki6_offset, Mki7_offset;
	    int firstcol = 0;
	    double M0[];
	    int M0_offset;
	    int k;

	    M0 = A;
	    M0_offset = A_offset;

	    while ( firstcol < n - 7 ) {        /* Do 8 columns */

	        Mki0 = M0;
	        Mki0_offset = M0_offset;
	        Mki1 = Mki0;
	        Mki1_offset = Mki0_offset + lda;
	        Mki2 = Mki1;
	        Mki2_offset = Mki1_offset + lda;
	        Mki3 = Mki2;
	        Mki3_offset = Mki2_offset + lda;
	        Mki4 = Mki3;
	        Mki4_offset = Mki3_offset + lda;
	        Mki5 = Mki4;
	        Mki5_offset = Mki4_offset + lda;
	        Mki6 = Mki5;
	        Mki6_offset = Mki5_offset + lda;
	        Mki7 = Mki6;
	        Mki7_offset = Mki6_offset + lda;

	        v00 = x0[x0_offset+firstcol];   v01 = x1[x1_offset+firstcol++];
	        v10 = x0[x0_offset+firstcol];   v11 = x1[x1_offset+firstcol++];
	        v20 = x0[x0_offset+firstcol];   v21 = x1[x1_offset+firstcol++];
	        v30 = x0[x0_offset+firstcol];   v31 = x1[x1_offset+firstcol++];
	        v40 = x0[x0_offset+firstcol];   v41 = x1[x1_offset+firstcol++];
	        v50 = x0[x0_offset+firstcol];   v51 = x1[x1_offset+firstcol++];
	        v60 = x0[x0_offset+firstcol];   v61 = x1[x1_offset+firstcol++];
	        v70 = x0[x0_offset+firstcol];   v71 = x1[x1_offset+firstcol++];

	        for (k = 0; k < m; k++) {
	            f0 = y0[y0_offset+k];
	            f1 = y1[y1_offset+k];
	            t0 = Mki0[Mki0_offset+k];  f0 += v00 * t0;  f1 += v01 * t0;
	            t1 = Mki1[Mki1_offset+k];  f0 += v10 * t1;  f1 += v11 * t1;
	            t2 = Mki2[Mki2_offset+k];  f0 += v20 * t2;  f1 += v21 * t2;
	            t3 = Mki3[Mki3_offset+k];  f0 += v30 * t3;  f1 += v31 * t3;
	            t4 = Mki4[Mki4_offset+k];  f0 += v40 * t4;  f1 += v41 * t4;
	            t5 = Mki5[Mki5_offset+k];  f0 += v50 * t5;  f1 += v51 * t5;
	            t6 = Mki6[Mki6_offset+k];  f0 += v60 * t6;  f1 += v61 * t6;
	            t7 = Mki7[Mki7_offset+k];  f0 += v70 * t7;  f1 += v71 * t7;
	            y0[y0_offset+k] = f0;
	            y1[y1_offset+k] = f1;
	        }

	        M0_offset += 8 * lda;
	    }

	    while ( firstcol < n - 3 ) {        /* Do 4 columns */
	        Mki0 = M0;
	        Mki0_offset = M0_offset;
	        Mki1 = Mki0;
	        Mki1_offset = Mki0_offset + lda;
	        Mki2 = Mki1;
	        Mki2_offset = Mki1_offset + lda;
	        Mki3 = Mki2;
	        Mki3_offset = Mki2_offset + lda;

	        v00 = x0[x0_offset+firstcol];   v01 = x1[x1_offset+firstcol++];
	        v10 = x0[x0_offset+firstcol];   v11 = x1[x1_offset+firstcol++];
	        v20 = x0[x0_offset+firstcol];   v21 = x1[x1_offset+firstcol++];
	        v30 = x0[x0_offset+firstcol];   v31 = x1[x1_offset+firstcol++];

	        for (k = 0; k < m; k++) {
	            f0 = y0[y0_offset+k];
	            f1 = y1[y1_offset+k];
	            t0 = Mki0[Mki0_offset+k];  f0 += v00 * t0;  f1 += v01 * t0;
	            t1 = Mki1[Mki1_offset+k];  f0 += v10 * t1;  f1 += v11 * t1;
	            t2 = Mki2[Mki2_offset+k];  f0 += v20 * t2;  f1 += v21 * t2;
	            t3 = Mki3[Mki3_offset+k];  f0 += v30 * t3;  f1 += v31 * t3;
	            y0[y0_offset+k] = f0;
	            y1[y1_offset+k] = f1;
	        }

	        M0_offset += 4 * lda;

	    }

	    while ( firstcol < n ) {            /* Do 1 column */
	        Mki0 = M0;
	        Mki0_offset = M0_offset;
	        v00 = x0[x0_offset+firstcol];   v01 = x1[x1_offset+firstcol++];

	        for (k = 0; k < m; k++) {
	            f0 = y0[y0_offset+k];
	            f1 = y1[y1_offset+k];
	            t0 = Mki0[Mki0_offset+k];
	            f0 += v00 * t0;
	            f1 += v01 * t0;
	            y0[y0_offset+k] = f0;
	            y1[y1_offset+k] = f1;
	        }

	        M0_offset += lda;
	    }

	}

}

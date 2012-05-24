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
package lu.jsuper;

public class Dlu_dmyblas2 {

	/**! \brief Solves a dense UNIT lower triangular system
	 *
	 *  The unit lower
	 * triangular matrix is stored in a 2D array M(1:nrow,1:ncol).
	 * The solution will be returned in the rhs vector.
	 */
	public static void dlsolve(int ldm, int ncol, double M[], double rhs[]) {
	    int k;
	    double x0, x1, x2, x3, x4, x5, x6, x7;
	    double[] M0;
	    int M0_offset;
	    double[] Mki0, Mki1, Mki2, Mki3, Mki4, Mki5, Mki6, Mki7;
	    int Mki0_offset, Mki1_offset, Mki2_offset, Mki3_offset, Mki4_offset,
	    Mki5_offset, Mki6_offset, Mki7_offset;
	    int firstcol = 0;

	    M0 = M;
	    M0_offset = 0;

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

		  x0 = rhs[firstcol];
		  x1 = rhs[firstcol+1] - x0 * Mki0[Mki0_offset];
		  Mki0_offset++;
		  x2 = rhs[firstcol+2] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset];
		  Mki0_offset++;
		  Mki1_offset++;
		  x3 = rhs[firstcol+3] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
				  - x2 * Mki2[Mki2_offset];
		  Mki0_offset++;
		  Mki1_offset++;
		  Mki2_offset++;
		  x4 = rhs[firstcol+4] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
				  - x2 * Mki2[Mki2_offset] - x3 * Mki3[Mki3_offset];
		  Mki0_offset++;
		  Mki1_offset++;
		  Mki2_offset++;
		  Mki3_offset++;
		  x5 = rhs[firstcol+5] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
				  - x2 * Mki2[Mki2_offset] - x3 * Mki3[Mki3_offset]
				  - x4 * Mki4[Mki4_offset];
		  Mki0_offset++;
		  Mki1_offset++;
		  Mki2_offset++;
		  Mki3_offset++;
		  Mki4_offset++;
		  x6 = rhs[firstcol+6] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
				  - x2 * Mki2[Mki2_offset] - x3 * Mki3[Mki3_offset]
				  - x4 * Mki4[Mki4_offset] - x5 * Mki5[Mki5_offset];
		  Mki0_offset++;
		  Mki1_offset++;
		  Mki2_offset++;
		  Mki3_offset++;
		  Mki4_offset++;
		  Mki5_offset++;
		  x7 = rhs[firstcol+7] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
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

		  rhs[++firstcol] = x1;
		  rhs[++firstcol] = x2;
		  rhs[++firstcol] = x3;
		  rhs[++firstcol] = x4;
		  rhs[++firstcol] = x5;
		  rhs[++firstcol] = x6;
		  rhs[++firstcol] = x7;
		  ++firstcol;

	      for (k = firstcol; k < ncol; k++) {
	    	  rhs[k] = rhs[k] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
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

	      x0 = rhs[firstcol];
	      x1 = rhs[firstcol+1] - x0 * Mki0[Mki0_offset];
	      Mki0_offset++;
	      x2 = rhs[firstcol+2] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset];
	      Mki0_offset++;
	      Mki1_offset++;
	      x3 = rhs[firstcol+3] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
	    		  - x2 * Mki2[Mki2_offset];
	      Mki0_offset++;
	      Mki1_offset++;
	      Mki2_offset++;

	      rhs[++firstcol] = x1;
	      rhs[++firstcol] = x2;
	      rhs[++firstcol] = x3;
	      ++firstcol;

	      for (k = firstcol; k < ncol; k++) {
	    	  rhs[k] = rhs[k] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset]
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

	      x0 = rhs[firstcol];
	      x1 = rhs[firstcol+1] - x0 * Mki0[Mki0_offset];
	      Mki0_offset++;

	      rhs[++firstcol] = x1;
	      ++firstcol;

	      for (k = firstcol; k < ncol; k++) {
	    	  rhs[k] = rhs[k] - x0 * Mki0[Mki0_offset] - x1 * Mki1[Mki1_offset];
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
	public static void dusolve(int ldm, int ncol, double M[], double rhs[])
	{
	    double xj;
	    int jcol, j, irow;

	    jcol = ncol - 1;

	    for (j = 0; j < ncol; j++) {

		xj = rhs[jcol] / M[jcol + jcol*ldm]; 		/* M(jcol, jcol) */
		rhs[jcol] = xj;

		for (irow = 0; irow < jcol; irow++)
		    rhs[irow] -= xj * M[irow + jcol*ldm];	/* M(irow, jcol) */

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
		double vec[],	/* in */
		double Mxvec[]	/* in/out */
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
	    M0_offset = 0;

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

		vi0 = vec[firstcol++];
		vi1 = vec[firstcol++];
		vi2 = vec[firstcol++];
		vi3 = vec[firstcol++];
		vi4 = vec[firstcol++];
		vi5 = vec[firstcol++];
		vi6 = vec[firstcol++];
		vi7 = vec[firstcol++];

		for (k = 0; k < nrow; k++) {
		    Mxvec[k] += vi0 * Mki0[Mki0_offset] + vi1 * Mki1[Mki1_offset]
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

		vi0 = vec[firstcol++];
		vi1 = vec[firstcol++];
		vi2 = vec[firstcol++];
		vi3 = vec[firstcol++];
		for (k = 0; k < nrow; k++) {
		    Mxvec[k] += vi0 * Mki0[Mki0_offset] + vi1 * Mki1[Mki1_offset]
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
		vi0 = vec[firstcol++];
		for (k = 0; k < nrow; k++) {
		    Mxvec[k] += vi0 * Mki0[Mki0_offset];
		    Mki0_offset++;
		}

		M0_offset += ldm;
	    }

	}

}

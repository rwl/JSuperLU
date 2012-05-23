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
/*
 * File name:		dmyblas2.c
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
	    double[] Mki0, Mki1, Mki2, Mki3, Mki4, Mki5, Mki6, Mki7;
	    int firstcol = 0;

	    M0 = &M[0];

	    while ( firstcol < ncol - 7 ) { /* Do 8 columns */
		  Mki0 = M0 + 1;
		  Mki1 = Mki0 + ldm + 1;
		  Mki2 = Mki1 + ldm + 1;
		  Mki3 = Mki2 + ldm + 1;
		  Mki4 = Mki3 + ldm + 1;
		  Mki5 = Mki4 + ldm + 1;
		  Mki6 = Mki5 + ldm + 1;
		  Mki7 = Mki6 + ldm + 1;

		  x0 = rhs[firstcol];
		  x1 = rhs[firstcol+1] - x0 * *Mki0++;
		  x2 = rhs[firstcol+2] - x0 * *Mki0++ - x1 * *Mki1++;
		  x3 = rhs[firstcol+3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++;
		  x4 = rhs[firstcol+4] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
		                   - x3 * *Mki3++;
		  x5 = rhs[firstcol+5] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
		                   - x3 * *Mki3++ - x4 * *Mki4++;
		  x6 = rhs[firstcol+6] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
		                   - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++;
		  x7 = rhs[firstcol+7] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
		                   - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++
				   - x6 * *Mki6++;

		  rhs[++firstcol] = x1;
		  rhs[++firstcol] = x2;
		  rhs[++firstcol] = x3;
		  rhs[++firstcol] = x4;
		  rhs[++firstcol] = x5;
		  rhs[++firstcol] = x6;
		  rhs[++firstcol] = x7;
		  ++firstcol;

	      for (k = firstcol; k < ncol; k++)
		rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
		                - x2 * *Mki2++ - x3 * *Mki3++
	                        - x4 * *Mki4++ - x5 * *Mki5++
				- x6 * *Mki6++ - x7 * *Mki7++;

	      M0 += 8 * ldm + 8;
	    }

	    while ( firstcol < ncol - 3 ) { /* Do 4 columns */
	      Mki0 = M0 + 1;
	      Mki1 = Mki0 + ldm + 1;
	      Mki2 = Mki1 + ldm + 1;
	      Mki3 = Mki2 + ldm + 1;

	      x0 = rhs[firstcol];
	      x1 = rhs[firstcol+1] - x0 * *Mki0++;
	      x2 = rhs[firstcol+2] - x0 * *Mki0++ - x1 * *Mki1++;
	      x3 = rhs[firstcol+3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++;

	      rhs[++firstcol] = x1;
	      rhs[++firstcol] = x2;
	      rhs[++firstcol] = x3;
	      ++firstcol;

	      for (k = firstcol; k < ncol; k++)
		rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
		                - x2 * *Mki2++ - x3 * *Mki3++;

	      M0 += 4 * ldm + 4;
	    }

	    if ( firstcol < ncol - 1 ) { /* Do 2 columns */
	      Mki0 = M0 + 1;
	      Mki1 = Mki0 + ldm + 1;

	      x0 = rhs[firstcol];
	      x1 = rhs[firstcol+1] - x0 * *Mki0++;

	      rhs[++firstcol] = x1;
	      ++firstcol;

	      for (k = firstcol; k < ncol; k++)
		rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++;

	    }
	}

}

/*! @file dpivotL.c
 * \brief Performs numerical pivoting
 *
 * <pre>
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 * Copyright (c) 1994 by Xerox Corporation.  All rights reserved.
 *
 * THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
 * EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
 *
 * Permission is hereby granted to use or copy this program for any
 * purpose, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is
 * granted, provided the above notices are retained, and a notice that
 * the code was modified is included with the above copyright notice.
 * </pre>
 */
package lu.jsuper;

import lu.jsuper.Dlu_slu_ddefs.GlobalLU_t;
import lu.jsuper.Dlu_slu_util.SuperLUStat_t;
import lu.jsuper.Dlu_superlu_enum_consts.PhaseType;

import static lu.jsuper.Dlu_util.DEBUG;
import static lu.jsuper.Dlu_util.MIN_COL;

import static lu.jsuper.Dlu_slu_util.EMPTY;


public class Dlu_dpivotL {

	/*! \brief
	 *
	 * <pre>
	 * Purpose
	 * =======
	 *   Performs the numerical pivoting on the current column of L,
	 *   and the CDIV operation.
	 *
	 *   Pivot policy:
	 *   (1) Compute thresh = u * max_(i>=j) abs(A_ij);
	 *   (2) IF user specifies pivot row k and abs(A_kj) >= thresh THEN
	 *           pivot row = k;
	 *       ELSE IF abs(A_jj) >= thresh THEN
	 *           pivot row = j;
	 *       ELSE
	 *           pivot row = m;
	 *
	 *   Note: If you absolutely want to use a given pivot order, then set u=0.0.
	 *
	 *   Return value: 0      success;
	 *                 i > 0  U(i,i) is exactly zero.
	 * </pre>
	 */
	@SuppressWarnings("unused")
	public static int dpivotL(
	        final int  jcol,      /* in */
	        final double u,       /* in - diagonal pivoting threshold */
	        int        usepr[],   /* re-use the pivot sequence given by perm_r/iperm_r */
	        int        perm_r[],  /* may be modified */
	        int        iperm_r[], /* in - inverse of perm_r */
	        int        iperm_c[], /* in - used to find diagonal of Pc*A*Pc' */
	        int        pivrow[],  /* out */
	        GlobalLU_t Glu,       /* modified - global LU data structures */
	        SuperLUStat_t stat    /* output */
			)
	{
	    int          fsupc;	    /* first column in the supernode */
	    int          nsupc;	    /* no of columns in the supernode */
	    int          nsupr;     /* no of rows in the supernode */
	    int          lptr;	    /* points to the starting subscript of the supernode */
	    int          pivptr, old_pivptr, diag, diagind;
	    double       pivmax, rtemp, thresh;
	    double       temp;
	    double       lu_sup_ptr[];
	    double       lu_col_ptr[];
	    int          lsub_ptr[];
	    int          isub, icol, k, itemp;
	    int          lsub[], xlsub[];
	    double       lusup[];
	    int          xlusup[];
	    float[]      ops = stat.ops;

	    /* Initialize pointers */
	    lsub       = Glu.lsub;
	    xlsub      = Glu.xlsub;
	    lusup      = Glu.lusup;
	    xlusup     = Glu.xlusup;
	    fsupc      = (Glu.xsup)[(Glu.supno)[jcol]];
	    nsupc      = jcol - fsupc;	        /* excluding jcol; nsupc >= 0 */
	    lptr       = xlsub[fsupc];
	    nsupr      = xlsub[fsupc+1] - lptr;
	    lu_sup_ptr = lusup;	/* start of the current supernode */
	    int lu_sup_ptr_offset = xlusup[fsupc];
	    lu_col_ptr = lusup;	/* start of jcol in the supernode */
	    int lu_col_ptr_offset = xlusup[jcol];
	    lsub_ptr   = lsub;	/* start of row indices of the supernode */
	    int lsub_ptr_offset = lptr;

		if (DEBUG) {
		if ( jcol == MIN_COL ) {
		    System.out.printf("Before cdiv: col %d\n", jcol);
		    for (k = nsupc; k < nsupr; k++)
		    	System.out.printf("  lu[%d] %f\n", lsub_ptr[lsub_ptr_offset + k],
		    			lu_col_ptr[lu_col_ptr_offset + k]);
		}
		}

	    /* Determine the largest abs numerical value for partial pivoting;
	       Also search for user-specified pivot, and diagonal element. */
	    if ( usepr[0] != 0 ) pivrow[0] = iperm_r[jcol];
	    diagind = iperm_c[jcol];
	    pivmax = 0.0;
	    pivptr = nsupc;
	    diag = EMPTY;
	    old_pivptr = nsupc;
	    for (isub = nsupc; isub < nsupr; ++isub) {
		rtemp = Math.abs (lu_col_ptr[lu_col_ptr_offset + isub]);
		if ( rtemp > pivmax ) {
		    pivmax = rtemp;
		    pivptr = isub;
		}
		if ( usepr[0] != 0 && lsub_ptr[lsub_ptr_offset + isub] == pivrow[0] )
			old_pivptr = isub;
		if ( lsub_ptr[lsub_ptr_offset + isub] == diagind ) diag = isub;
	    }

	    /* Test for singularity */
	    if ( pivmax == 0.0 ) {
	    if (true) {
			pivrow[0] = lsub_ptr[lsub_ptr_offset + pivptr];
			perm_r[pivrow[0]] = jcol;
	    } else {
	    	perm_r[diagind] = jcol;
	    }
		usepr[0] = 0;
		return (jcol+1);
	    }

	    thresh = u * pivmax;

	    /* Choose appropriate pivotal element by our policy. */
	    if ( usepr[0] != 0 ) {
	        rtemp = Math.abs (lu_col_ptr[lu_col_ptr_offset + old_pivptr]);
		if ( rtemp != 0.0 && rtemp >= thresh )
		    pivptr = old_pivptr;
		else
		    usepr[0] = 0;
	    }
	    if ( usepr[0] == 0 ) {
		/* Use diagonal pivot? */
		if ( diag >= 0 ) { /* diagonal exists */
		    rtemp = Math.abs (lu_col_ptr[lu_col_ptr_offset + diag]);
		    if ( rtemp != 0.0 && rtemp >= thresh ) pivptr = diag;
	        }
		pivrow[0] = lsub_ptr[lsub_ptr_offset + pivptr];
	    }

	    /* Record pivot row */
	    perm_r[pivrow[0]] = jcol;

	    /* Interchange row subscripts */
	    if ( pivptr != nsupc ) {
		itemp = lsub_ptr[lsub_ptr_offset + pivptr];
		lsub_ptr[lsub_ptr_offset + pivptr] = lsub_ptr[lsub_ptr_offset + nsupc];
		lsub_ptr[lsub_ptr_offset + nsupc] = itemp;

		/* Interchange numerical values as well, for the whole snode, such
		 * that L is indexed the same way as A.
	 	 */
		for (icol = 0; icol <= nsupc; icol++) {
		    itemp = pivptr + icol * nsupr;
		    temp = lu_sup_ptr[lu_sup_ptr_offset + itemp];
		    lu_sup_ptr[lu_sup_ptr_offset + itemp] = lu_sup_ptr[lu_sup_ptr_offset + nsupc + icol*nsupr];
		    lu_sup_ptr[lu_sup_ptr_offset + nsupc + icol*nsupr] = temp;
		}
	    } /* if */

	    /* cdiv operation */
	    ops[PhaseType.FACT.ordinal()] += nsupr - nsupc;

	    temp = 1.0 / lu_col_ptr[lu_col_ptr_offset + nsupc];
	    for (k = nsupc+1; k < nsupr; k++)
		lu_col_ptr[lu_col_ptr_offset + k] *= temp;

	    return 0;
	}

}

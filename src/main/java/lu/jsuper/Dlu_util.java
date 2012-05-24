/*! @file Dlu_util.java
 * \brief Utility functions
 *
 * <pre>
 * -- SuperLU routine (version 4.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November, 2010
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

import lu.jsuper.Dlu_slu_util.Dlu_superlu_options_t;
import lu.jsuper.Dlu_slu_util.SuperLUStat_t;
import lu.jsuper.Dlu_superlu_enum_consts.IterRefine_t;
import lu.jsuper.Dlu_superlu_enum_consts.PhaseType;
import lu.jsuper.Dlu_superlu_enum_consts.colperm_t;
import lu.jsuper.Dlu_superlu_enum_consts.fact_t;
import lu.jsuper.Dlu_superlu_enum_consts.trans_t;
import lu.jsuper.Dlu_superlu_enum_consts.yes_no_t;

import static lu.jsuper.Dlu_sp_ienv.sp_ienv;
import static lu.jsuper.Dlu_slu_util.SUPERLU_MAX;
import static lu.jsuper.Dlu_memory.intCalloc;


public class Dlu_util {

	public static int PRNTlevel = 0;
	public static boolean DEBUG = false;
	public static boolean USE_VENDOR_BLAS = false;
	public static int MIN_COL = 8;

	public static void superlu_abort_and_exit(String msg) {
		System.err.print(msg);
		System.exit(-1);
	}

	/*! \brief Set the default values for the options argument.
	 */
	public static void set_default_options(Dlu_superlu_options_t options) {
	    options.Fact = fact_t.DOFACT;
	    options.Equil = yes_no_t.YES;
	    options.ColPerm = colperm_t.COLAMD;
	    options.Trans = trans_t.NOTRANS;
	    options.IterRefine = IterRefine_t.NOREFINE;
	    options.DiagPivotThresh = 1.0;
	    options.SymmetricMode = yes_no_t.NO;
	    options.PivotGrowth = yes_no_t.NO;
	    options.ConditionNumber = yes_no_t.NO;
	    options.PrintStat = yes_no_t.YES;
	}

	public static void StatInit(SuperLUStat_t stat) {
		int i, w, panel_size, relax;

	    panel_size = sp_ienv(1);
	    relax = sp_ienv(2);
	    w = SUPERLU_MAX(panel_size, relax);
	    stat.panel_histo = intCalloc(w+1);
	    stat.utime = new double[PhaseType.NPHASES.ordinal()];
	    stat.ops = new float[PhaseType.NPHASES.ordinal()];
	    for (i = 0; i < PhaseType.NPHASES.ordinal(); ++i) {
	        stat.utime[i] = 0.;
	        stat.ops[i] = 0.f;
	    }
	    stat.TinyPivots = 0;
	    stat.RefineSteps = 0;
	    stat.expansions = 0;
		if ( PRNTlevel >= 1 ) {
		    System.out.printf(".. parameters in sp_ienv():\n");
		    System.out.printf("\t 1: panel size \t %4d \n" +
		           "\t 2: relax      \t %4d \n" +
		           "\t 3: max. super \t %4d \n" +
		           "\t 4: row-dim 2D \t %4d \n" +
		           "\t 5: col-dim 2D \t %4d \n" +
		           "\t 6: fill ratio \t %4d \n",
			   sp_ienv(1), sp_ienv(2), sp_ienv(3),
			   sp_ienv(4), sp_ienv(5), sp_ienv(6));
		}
	}

	/**! \brief Fills an integer array with a given value.
	 */
	public static void ifill(int a[], int alen, int ival)
	{
	    int i;
	    for (i = 0; i < alen; i++) a[i] = ival;
	}

}

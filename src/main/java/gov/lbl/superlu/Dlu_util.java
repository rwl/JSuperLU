/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 */
package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_pdsp_defs.pxgstrf_shared_t;
import gov.lbl.superlu.Dlu_slu_mt_util.Gstat_t;
import gov.lbl.superlu.Dlu_slu_mt_util.cp_panel_t;
import gov.lbl.superlu.Dlu_slu_mt_util.desc_eft_t;
import gov.lbl.superlu.Dlu_slu_mt_util.panstat_t;
import gov.lbl.superlu.Dlu_slu_mt_util.procstat_t;
import gov.lbl.superlu.Dlu_supermatrix.DNformat;
import gov.lbl.superlu.Dlu_supermatrix.NCPformat;
import gov.lbl.superlu.Dlu_supermatrix.NCformat;
import gov.lbl.superlu.Dlu_supermatrix.SCPformat;
import gov.lbl.superlu.Dlu_supermatrix.SCformat;
import gov.lbl.superlu.Dlu_supermatrix.SuperMatrix;

import static gov.lbl.superlu.Dlu.DEBUGlevel;
import static gov.lbl.superlu.Dlu.GEMV2;
import static gov.lbl.superlu.Dlu.PREDICT_OPT;
import static gov.lbl.superlu.Dlu.PRNTlevel;
import static gov.lbl.superlu.Dlu.PROFILE;
import static gov.lbl.superlu.Dlu.SCATTER_FOUND;
import static gov.lbl.superlu.Dlu.USE_VENDOR_BLAS;
import static gov.lbl.superlu.Dlu.exit;
import static gov.lbl.superlu.Dlu.fflush;
import static gov.lbl.superlu.Dlu.fprintf;
import static gov.lbl.superlu.Dlu.printf;
import static gov.lbl.superlu.Dlu.stderr;
import static gov.lbl.superlu.Dlu.stdout;
import static gov.lbl.superlu.Dlu_pxgstrf_synch.panel_t.TREE_DOMAIN;
import static gov.lbl.superlu.Dlu_slu_mt_util.EMPTY;
import static gov.lbl.superlu.Dlu_slu_mt_util.SINGLETON;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_ABORT;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_MAX;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPER_REP;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.FACT;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.NPHASES;
import static gov.lbl.superlu.Dlu_slu_mt_util.PhaseType.SOLVE;
import static gov.lbl.superlu.Dlu_slu_mt_util.how_selected_t.DADPAN;
import static gov.lbl.superlu.Dlu_slu_mt_util.how_selected_t.NOPIPE;
import static gov.lbl.superlu.Dlu_slu_mt_util.how_selected_t.PIPE;




public class Dlu_util {

	static
	void superlu_abort_and_exit(String msg)
	{
	    fprintf(stderr, msg);
	    exit (-1);
	}

	/* Deallocate the structure pointing to the actual storage of the matrix. */
	static
	void
	Destroy_SuperMatrix_Store(SuperMatrix A)
	{
	    A.Store = null;
	if ( PRNTlevel==1 ) {
	    printf(".. Destroy_SuperMatrix_Store ...\n");
	}
	}

	/* A is of type Stype==NC */
	static
	void
	Destroy_CompCol_Matrix(SuperMatrix A)
	{
	    ((NCformat)A.Store).rowind = null;
	    ((NCformat)A.Store).colptr = null;
	    ((NCformat)A.Store).nzval = null;
	    A.Store = null;
	if ( PRNTlevel==1 ) {
	    printf(".. Destroy_CompCol_Matrix ...\n");
	}
	}

	/* A is of type Stype==NCP */
	static
	void
	Destroy_CompCol_Permuted(SuperMatrix A)
	{
	    ((NCPformat)A.Store).colbeg = null;
	    ((NCPformat)A.Store).colend = null;
	    A.Store = null;
	if ( PRNTlevel==1 ) {
	    printf(".. Destroy_CompCol_Permuted ...\n");
	}
	}

	/* A is of type Stype==NCP */
	static
	void
	Destroy_CompCol_NCP(SuperMatrix A)
	{
	    ((NCPformat)A.Store).nzval = null;
	    ((NCPformat)A.Store).rowind = null;
	    ((NCPformat)A.Store).colbeg = null;
	    ((NCPformat)A.Store).colend = null;
	    A.Store = null;
	if ( PRNTlevel==1 ) {
	    printf(".. Destroy_CompCol_NCP ...\n");
	}
	}

	/* A is of type Stype==SC */
	static
	void
	Destroy_SuperNode_Matrix(SuperMatrix A)
	{
	    ((SCformat)A.Store).rowind = null;
	    ((SCformat)A.Store).rowind_colptr = null;
	    ((SCformat)A.Store).nzval = null;
	    ((SCformat)A.Store).nzval_colptr = null;
	    ((SCformat)A.Store).col_to_sup = null;
	    ((SCformat)A.Store).sup_to_col = null;
	    A.Store = null;
	if ( PRNTlevel==1 ) {
	    printf(".. Destroy_SuperNode_Matrix ...\n");
	}
	}

	/* A is of type Stype==SCP */
	static
	void
	Destroy_SuperNode_SCP(SuperMatrix A)
	{
	    ((SCPformat)A.Store).rowind = null;
	    ((SCPformat)A.Store).rowind_colbeg = null;
	    ((SCPformat)A.Store).rowind_colend = null;
	    ((SCPformat)A.Store).nzval = null;
	    ((SCPformat)A.Store).nzval_colbeg = null;
	    ((SCPformat)A.Store).nzval_colend = null;
	    ((SCPformat)A.Store).col_to_sup = null;
	    ((SCPformat)A.Store).sup_to_colbeg = null;
	    ((SCPformat)A.Store).sup_to_colend = null;
	    A.Store = null;
	if ( PRNTlevel==1 ) {
	    printf(".. Destroy_SuperNode_SCP ...\n");
	}
	}

	/*
	 * Reset repfnz[*] for the current column
	 */
	static
	void
	pxgstrf_resetrep_col(final int nseg, final int segrep[], int repfnz[])
	{
	    int i, irep;

	    for (i = 0; i < nseg; ++i) {
		irep = segrep[i];
		repfnz[irep] = EMPTY;
	    }
	}


	/*
	 * Count the total number of nonzeros in factors L and U,  and in the
	 * symmetrically reduced L.
	 */
	static
	void
	countnz(final int n, int xprune[], int nnzL[], int nnzU[], GlobalLU_t Glu)
	{
	    int nsuper, fsupc, i, j, nnzL0, jlen, irep;
	    int nnzsup = 0;
	    int[] xsup, xsup_end, xlsub, xlsub_end, supno;

	    xsup      = Glu.xsup;
	    xsup_end  = Glu.xsup_end;
	    xlsub     = Glu.xlsub;
	    xlsub_end = Glu.xlsub_end;
	    supno     = Glu.supno;
	    nnzU[0]   = Glu.nextu;
	    nnzL0     = 0;
	    nnzL[0]   = 0;
	    nsuper    = supno[n];

	    if ( n <= 0 ) return;

	    /*
	     * For each supernode ...
	     */
	    for (i = 0; i <= nsuper; i++) {
		fsupc = xsup[i];
		jlen = xlsub_end[fsupc] - xlsub[fsupc];
		nnzsup += jlen * (xsup_end[i] - fsupc);

		for (j = fsupc; j < xsup_end[i]; j++) {
		    nnzL[0] += jlen;
		    nnzU[0] += j - fsupc + 1;
		    jlen--;
		}
		irep = SUPER_REP(xsup_end, i);
		if ( SINGLETON(xsup_end, xsup_end, supno[irep]) )
		    nnzL0 += xprune[irep] - xlsub_end[irep];
		else
		    nnzL0 += xprune[irep] - xlsub[irep];
	    }

	if ( PRNTlevel==1 ) {
	    printf(".. # supernodes = %d\n", nsuper+1);
	    printf(".. # edges in symm-reduced L = %d\n", nnzL0);
	    if ( Glu.dynamic_snode_bound != 0 )
	      printf(".. # NZ in LUSUP %d, dynamic bound %d, utilization %.2f\n",
		     nnzsup, Glu.nextlu, (float)nnzsup/Glu.nextlu);
	    else
	      printf(".. # NNZ in LUSUP %d, static bound %d, utilization %.2f\n",
		     nnzsup, Glu.nzlumax, (float)nnzsup/Glu.nzlumax);
	}
	}



	/*
	 * Fix up the data storage lsub for L-subscripts. It reclaims the
	 * storage for the adjancency lists of the pruned graph, and applies
	 * row permuation to the row subscripts of matrix $L$.
	 */
	static
	void
	fixupL(final int n, final int perm_r[], GlobalLU_t Glu)
	{
	    int nsuper, fsupc, nextl, i, j, jstrt;
	    int[] xsup, xsup_end, lsub, xlsub, xlsub_end;

	    if ( n <= 1 ) return;

	    xsup      = Glu.xsup;
	    xsup_end  = Glu.xsup_end;
	    lsub      = Glu.lsub;
	    xlsub     = Glu.xlsub;
	    xlsub_end = Glu.xlsub_end;
	    nsuper    = Glu.supno[n];
	    nextl     = 0;

	    /*
	     * For each supernode ...
	     */
	    for (i = 0; i <= nsuper; i++) {
		fsupc = xsup[i];
		jstrt = xlsub[fsupc];
		xlsub[fsupc] = nextl;
		for (j = jstrt; j < xlsub_end[fsupc]; j++) {
		    lsub[nextl] = perm_r[lsub[j]]; /* Now indexed into P*A */
		    nextl++;
	  	}
		xlsub_end[fsupc] = nextl;
	    }
	    xlsub[n] = nextl;

	if ( PRNTlevel==1 ) {
	    printf(".. # edges in supernodal graph of L = %d\n", nextl);
	    fflush(stdout);
	}
	}

	/*
	 * Print all definitions to be used by CPP.
	 */
	static
	int cpp_defs()
	{
	    printf("CPP Defs:\n");
	if (true) {
	    printf("\tPRNTlevel=%d\n", PRNTlevel);
	}
	if (true) {
	    printf("\tDEBUGlevel=%d\n", DEBUGlevel);
	}
	if (PROFILE) {
	    printf("\tPROFILE\n");
	}
	if (PREDICT_OPT) {
	    printf("\tPREDICT_OPT\n");
	}
	if (USE_VENDOR_BLAS) {
	    printf("\tUSE_VENDOR_BLAS\n");
	}
	if (GEMV2) {
	    printf("\tGEMV2\n");
	}
	if (SCATTER_FOUND) {
	    printf("\tSCATTER_FOUND\n");
	}

	    return 0;
	}

	/*
	 * Compress the data storage LUSUP[*] for supernodes. It removes the
	 * memory holes due to untightness of the upper bounds by A'A-supernode.
	 */
	static
	void
	compressSUP(final int n, GlobalLU_t Glu)
	{
	    int nextlu, i, j, jstrt;
	    int xlusup[], xlusup_end[];
	    double lusup[];

	    if ( n <= 1 ) return;

	    lusup     = Glu.lusup;
	    xlusup    = Glu.xlusup;
	    xlusup_end= Glu.xlusup_end;
	    nextlu     = 0;

	    for (j = 0; j < n; ++j) {
		jstrt = xlusup[j];
		xlusup[j] = nextlu;
		for (i = jstrt; i < xlusup_end[j]; ++i, ++nextlu)
		    lusup[nextlu] = lusup[i];
		xlusup_end[j] = nextlu;
	    }
	    xlusup[n] = nextlu;
	    printf("\tcompressSUP() nextlu %d\n", nextlu);
	}

	static
	int check_mem_leak(String where)
	{
//	    void *addr;
//	    addr = (void *)sbrk(0);
//	    printf("\tsbrk(0) %s: addr = %u\n", where, addr);
	    return 0;
	}

	/*
	 * Diagnostic print of segment info after pdgstrf_panel_dfs().
	 */
	static
	void print_panel_seg(int n, int w, int jcol, int nseg,
			     int segrep[], int repfnz[])
	{
	    int j, k;

	    for (j = jcol; j < jcol+w; j++) {
		printf("\tcol %d:\n", j);
		for (k = 0; k < nseg; k++)
		    printf("\t\tseg %d, segrep %d, repfnz %d\n", k,
				segrep[k], repfnz[(j-jcol)*n + segrep[k]]);
	    }
	}

	/*
	 * Allocate storage for various statistics.
	 */
	static
	void
	StatAlloc(final int n, final int nprocs, final int panel_size,
			final int relax, Gstat_t Gstat)
	{
	    int w;

	    w = SUPERLU_MAX( panel_size, relax ) + 1;
	    Gstat.panel_histo = intCalloc(w);
	    Gstat.utime = (double []) doubleMalloc(NPHASES);
	    Gstat.ops   = new float[NPHASES.ordinal()];

	    if ( (Gstat.procstat = new procstat_t[nprocs]) == null )
		SUPERLU_ABORT( "SUPERLU_MALLOC failed for procstat[]" );

	if (PRNTlevel==1) {
	    printf(".. StatAlloc(): n %d, nprocs %d, panel_size %d, relax %d\n",
			n, nprocs, panel_size, relax);
	}
	if (PROFILE) {
	    if ( (Gstat.panstat = new panstat_t[n]) == null )
		SUPERLU_ABORT( "SUPERLU_MALLOC failed for panstat[]" );
	    Gstat.panhows = intCalloc(3);
	    Gstat.height = intCalloc(n+1);
	    if ( (Gstat.flops_by_height = new float[n]) == null )
		SUPERLU_ABORT("SUPERLU_MALLOC failed for flops_by_height[]");
	}

	if (PREDICT_OPT) {
	    if ( (Gstat.cp_panel = new cp_panel_t[n]) == null )
		SUPERLU_ABORT( "SUPERLU_MALLOC failed for cp_panel[]" );
	    if ( (Gstat.desc_eft = new desc_eft_t[n]) == null )
		SUPERLU_ABORT( "SUPERLU_MALLOC failed for desc_eft[]" );
	    Gstat.cp_firstkid = intMalloc(n+1);
	    Gstat.cp_nextkid = intMalloc(n+1);
	}

	}

	/*
	 * Initialize various statistics variables.
	 */
	static
	void
	StatInit(final int n, final int nprocs, Gstat_t Gstat)
	{
	    int i;

	    for (i = 0; i < NPHASES.ordinal(); ++i) {
		Gstat.utime[i] = 0;
		Gstat.ops[i] = 0;
	    }

	    for (i = 0; i < nprocs; ++i) {
		Gstat.procstat[i].panels = 0;
		Gstat.procstat[i].fcops = 0.0f;
		Gstat.procstat[i].skedwaits = 0;
		Gstat.procstat[i].skedtime = 0.0;
		Gstat.procstat[i].cs_time = 0.0;
		Gstat.procstat[i].spintime = 0.0;
		Gstat.procstat[i].pruned = 0;
		Gstat.procstat[i].unpruned = 0;
	    }

	if (PROFILE) {
	    for (i = 0; i < n; ++i) {
		Gstat.panstat[i].fctime = 0.0;
		Gstat.panstat[i].flopcnt = 0.0f;
		Gstat.panstat[i].pipewaits = 0;
		Gstat.panstat[i].spintime = 0.0;
		Gstat.flops_by_height[i] = 0.0f;
	    }
	    for (i = 0; i < 3; ++i) Gstat.panhows[i] = 0;
	    Gstat.dom_flopcnt = 0.f;
	    Gstat.flops_last_P_panels = 0;
	}

	if (PREDICT_OPT) {
	    for (i = 0; i < n; ++i)
	    Gstat.cp_panel[i].est = Gstat.cp_panel[i].pdiv = 0;
	}

	if ( PRNTlevel==1 ) {
	    printf(".. StatInit(): n %d, nprocs %d\n", n, nprocs);
	}
	}


	/* Print timings used in factorization and solve. */
	static
	void
	PrintStat(Gstat_t Gstat)
	{
	    double         utime[];
	    float        ops[];

	    utime = Gstat.utime;
	    ops   = Gstat.ops;
	    printf("Factor time  = %8.2f\n", utime[FACT.ordinal()]);
	    if ( utime[FACT.ordinal()] != 0.0 )
	      printf("Factor flops = %e\tMflops = %8.2f\n", ops[FACT.ordinal()],
		     ops[FACT.ordinal()]*1e-6/utime[FACT.ordinal()]);

	    printf("Solve time   = %8.2f\n", utime[SOLVE.ordinal()]);
	    if ( utime[SOLVE.ordinal()] != 0.0 )
	      printf("Solve flops = %e\tMflops = %8.2f\n", ops[SOLVE.ordinal()],
		     ops[SOLVE.ordinal()]*1e-6/utime[SOLVE.ordinal()]);

	}

	static
	void
	StatFree(Gstat_t Gstat)
	{
	    Gstat.panel_histo = null;
	    Gstat.utime = null;
	    Gstat.ops = null;
	    Gstat.procstat = null;

	if (PROFILE) {
	    Gstat.panstat = null;
	    Gstat.panhows = null;
	    Gstat.height = null;
	    Gstat.flops_by_height = null;
	}

	if (PREDICT_OPT) {
	    Gstat.cp_panel = null;
	    Gstat.desc_eft = null;
	    Gstat.cp_firstkid = null;
	    Gstat.cp_nextkid = null;
	}

	if (PRNTlevel==1) {
	    printf(".. StatFree(): Free Stat variables.\n");
	}
	}

	static
	float
	LUFactFlops(Gstat_t Gstat)
	{
	    return (Gstat.ops[FACT.ordinal()]);
	}

	static
	float
	LUSolveFlops(Gstat_t Gstat)
	{
	    return (Gstat.ops[SOLVE.ordinal()]);
	}



	/*
	 * Fills an integer array with a given value.
	 */
	static
	void ifill(int a[], int alen, int ival)
	{
	    int i;
	    for (i = 0; i < alen; i++) a[i] = ival;
	}



	/*
	 * Get the statistics of the supernodes
	 */
	static int NBUCKS = 10;
	static int	max_sup_size;

	static
	void super_stats(int nsuper, int xsup[], int xsup_end[])
	{
	    int nsup1 = 0;
	    int          i, isize, whichb, bl, bh;
	    int          bucket[] = new int[NBUCKS];

	    max_sup_size = 0;

	    /* Histogram of the supernode sizes */
	    ifill (bucket, NBUCKS, 0);

	    for (i = 0; i <= nsuper; i++) {
	        isize = xsup_end[i] - xsup[i];
		if ( isize == 1 ) nsup1++;
		if ( max_sup_size < isize ) max_sup_size = isize;
	        whichb = (int) (float) isize / max_sup_size * NBUCKS;
	        if (whichb >= NBUCKS) whichb = NBUCKS - 1;
	        bucket[whichb]++;
	    }

	    printf("** Supernode statistics:\n\tno of supernodes = %d\n", nsuper+1);
	    printf("\tmax supernode size = %d\n", max_sup_size);
	    printf("\tno of size 1 supernodes = %d\n", nsup1);

	    printf("\tHistogram of supernode size:\n");
	    for (i = 0; i < NBUCKS; i++) {
	        bl = (int) (float) i * max_sup_size / NBUCKS;
	        bh = (int) (float) (i+1) * max_sup_size / NBUCKS;
	        printf("\t%3d-%3d\t\t%d\n", bl+1, bh, bucket[i]);
	    }

	}

	static
	void panel_stats(int n, int max_w, int in_domain[], Gstat_t Gstat)
	{
	    int i, w;
	    float histo_flops[], total;

	    histo_flops = new float[max_w];

	    for (i = 0; i < max_w; ++i) histo_flops[i] = 0;
	    total = 0;
	    for (i = 0; i < n; i += w) {
		w = Gstat.panstat[i].size;
		if ( in_domain[i] != TREE_DOMAIN.ordinal() ) {
		    histo_flops[w - 1] += Gstat.panstat[i].flopcnt;
		    total += Gstat.panstat[i].flopcnt;
		}
	    }

	    if ( total != 0.0 ) {
		printf("** Panel & flops distribution: nondomain flopcnt %e\n", total);
		for (i = 1; i <= max_w; i++)
		    printf("\t%d\t%d\t%e (%.2f)\n", i, Gstat.panel_histo[i],
			   histo_flops[i-1], histo_flops[i-1]/total);
	    }
	}



	static
	float SpaSize(int n, int np, float sum_npw)
	{
	    return (sum_npw*8 + np*8 + n*4)/1024.f;
	}

	static
	float DenseSize(int n, float sum_nw)
	{
	    return (sum_nw*8 + n*8)/1024.f;
	}


	/*
	 * Check whether repfnz[] == EMPTY after reset.
	 */
	void check_repfnz(int n, int w, int jcol, int repfnz[])
	{
	    int jj, k;

	    for (jj = jcol; jj < jcol+w; jj++)
		for (k = 0; k < n; k++)
		    if ( repfnz[(jj-jcol)*n + k] != EMPTY ) {
			fprintf(stderr, "col %d, repfnz_col[%d] = %d\n", jj,
				k, repfnz[(jj-jcol)*n + k]);
			SUPERLU_ABORT("repfnz[] not empty.");
		    }
	}


	static
	int PrintInt10(String name, int len, int x[])
	{
	    int i;

	    printf("(len=%d) %s:", len, name);
	    for (i = 0; i < len; ++i) {
		if ( i % 10 == 0 ) printf("\n[%4d-%4d]", i, i+9);
		printf("%6d", x[i]);
	    }
	    printf("\n");
	    return 0;
	}

	/* Print a summary of the testing results. */
	static
	void
	PrintSumm(String type, int nfail, int nrun, int nerrs)
	{
	    if ( nfail > 0 )
		printf("%3s driver: %d out of %d tests failed to pass the threshold\n",
		       type, nfail, nrun);
	    else
		printf("All tests for %3s driver passed the threshold (%6d tests run)\n", type, nrun);

	    if ( nerrs > 0 )
		printf("%6d error messages recorded\n", nerrs);
	}


	/* Print the adjacency list for graph of L, including the pruned graph,
	   graph of U, and L supernodes partition */
	static
	int PrintGLGU(int n, int xprune[], GlobalLU_t Glu)
	{
	    int nsuper = Glu.nsuper;
	    PrintInt10("LSUB", Glu.xlsub_end[n-1], Glu.lsub);
	    PrintInt10("XLSUB", n, Glu.xlsub);
	    PrintInt10("XLSUB_END", n, Glu.xlsub_end);
	    PrintInt10("XPRUNE", n, xprune);
	    PrintInt10("USUB", Glu.xusub_end[n-1], Glu.usub);
	    PrintInt10("XUSUB", n, Glu.xusub);
	    PrintInt10("XUSUB_END", n, Glu.xusub_end);
	    PrintInt10("SUPNO", n, Glu.supno);
	    PrintInt10("XSUP", nsuper+1, Glu.xsup);
	    PrintInt10("XSUP_END", nsuper+1, Glu.xsup_end);
	    return 0;
	}

//	/*
//	 * Print the statistics of the relaxed snodes for matlab process
//	 */
//	static
//	void relax_stats(int start, int end, int step)
//	{
//	    FILE fp;
//	    int i;
//
//	    fp = fopen("relax.m", "w");
//
//	    fprintf(fp,"relax = [\n");
//	    for (i = start; i <= end; i += step) fprintf(fp, "%d ", i);
//	    fprintf(fp, "];\n");
//
//	    fprintf(fp, "fctime = [\n");
//	    for (i = start; i <= end; i += step)
//		fprintf(fp, "%15.8e\n ", stat_relax[i].fctime);
//	    fprintf(fp, "];\n");
//
//	    fprintf(fp, "mflops = [\n");
//	    for (i = start; i <= end; i += step)
//		fprintf(fp, "%15.8e\n ", (float)stat_relax[i].flops / 1e6);
//	    fprintf(fp, "];\n");
//
//	    fprintf(fp, "mnzs = [\n");
//	    for (i = start; i <= end; i += step)
//	 	fprintf(fp, "%15.8e\n ", stat_relax[i].nzs / 1e6);
//	    fprintf(fp, "];\n");
//
//	    fclose(fp);
//	}
//
//	/*
//	 * Obtain the distribution of time/flops/nzs on the snode size.
//	 */
//	static
//	void snode_profile(int nsuper, int xsup[])
//	{
//	    FILE fp;
//	    int i, j;
//	    int ssize;
//
//	    if ( !(stat_snode = (stat_snode_t *) SUPERLU_MALLOC((max_sup_size+1) *
//		sizeof(stat_snode_t))) ) ABORT("SUPERLU_MALLOC fails for stat_snode[].");
//
//	    for (i = 0; i <= max_sup_size; i++) {
//		stat_snode[i].ncols = 0;
//		stat_snode[i].flops = 0;
//		stat_snode[i].nzs = 0;
//		stat_snode[i].fctime = 0.0;
//	    }
//
//	    for (i = 0; i <= nsuper; i++) {
//
//		ssize = xsup[i+1] - xsup[i];
//		stat_snode[ssize].ncols += ssize;
//
//	        for (j=xsup[i]; j<xsup[i+1]; j++) {
//		    stat_snode[ssize].flops += stat_col[j].flops;
//		    stat_snode[ssize].nzs += stat_col[j].nzs;
//		    stat_snode[ssize].fctime += stat_col[j].fctime;
//		}
//
//	    }
//
//	    fp = fopen("snode.m", "w");
//
//	    fprintf(fp, "max_sup_size = %d;\n", max_sup_size);
//
//	    fprintf(fp,"ncols = [");
//	    for (i = 1; i <= max_sup_size; i++)
//		fprintf(fp, "%d ", stat_snode[i].ncols);
//	    fprintf(fp, "];\n");
//
//	    fprintf(fp, "fctime = [");
//	    for (i = 1; i <= max_sup_size; i++)
//		fprintf(fp, "%15.8e\n", stat_snode[i].fctime);
//	    fprintf(fp, "];\n");
//
//	    fprintf(fp, "mflops = [");
//	    for (i = 1; i <= max_sup_size; i++)
//		fprintf(fp, "%15.8e\n", (float) stat_snode[i].flops / 1e6);
//	    fprintf(fp, "];\n");
//
//	    fprintf(fp, "mnzs = [");
//	    for (i = 1; i <= max_sup_size; i++)
//		fprintf(fp, "%15.8e\n", (float) stat_snode[i].nzs / 1e6);
//	    fprintf(fp, "];\n");
//
//	    fclose(fp);
//
//	    SUPERLU_FREE (stat_snode);
//
//	}

	static
	int print_int_vec(String what, int n, int vec[])
	{
	    int i;
	    printf("%s\n", what);
	    for (i = 0; i < n; ++i) printf("%d\t%d\n", i, vec[i]);
	    return 0;
	}


	/*
	 * Print the parallel execution statistics.
	 */
	static
	int ParallelProfile(final int n, final int supers, final int panels,
			final int procs, Gstat_t Gstat)
	{
	    int i, imax, pruned, unpruned, waits, itemp, cs_numbers;
	    float loadmax, loadtot, temp, thresh, loadprint;
	    float waittime, cs_time;
	    double    utime[] = Gstat.utime;
	    procstat_t pstat;
	    panstat_t pan;
//	    void print_flops_by_height(int, panstat_t *, int *, float *);

	    printf("\n---- Parallel Profile Per Processor ----\n");
	    printf("%4s%16s%8s%10s%10s%10s%10s%8s\n", "proc", "factops",
		   "seconds", "skedwaits", "skedtime", "CS-time",
		   "spin-time", "[%tot]");
	    for (i = 0; i < procs; ++i) {
		pstat = (Gstat.procstat[i]);
		if ( pstat.fctime != 0 ) {
		    temp = (float) (pstat.spintime/pstat.fctime*100.f);
		    printf("%4d%16e%8.2f%10d%10.3f%10.3f%10.3f%8.1f\n",
			   i, pstat.fcops, pstat.fctime, pstat.skedwaits,
			   pstat.skedtime, pstat.cs_time, pstat.spintime, temp);
		}
	    }

	    printf("%4s%8s%12s%14s\n",
		   "proc", "#panels", "dfs_pruned","dfs_unpruned");
	    pruned = unpruned = 0;
	    cs_time = 0.0f;
	    for (i = 0; i < procs; ++i) {
		pstat = (Gstat.procstat[i]);
		printf("%4d%8d%12d%14d\n", i, pstat.panels,
			pstat.pruned, pstat.unpruned);
		pruned += Gstat.procstat[i].pruned;
		unpruned += Gstat.procstat[i].unpruned;
		cs_time += Gstat.procstat[i].cs_time;
	    }
	    temp = pruned + unpruned;
	    if ( temp != 0 ) {
	    	printf("%12s%26s\n", "", "--------------------");
	    	printf("%12s%12d%14d%14.0f\n", "total", pruned, unpruned, temp);
	    	printf("%12s%12.2f%14.2f\n", "frac.", pruned/temp, unpruned/temp);
	    }

	    printf("%16s%16d\n", "piped-panels", Gstat.panhows[PIPE.ordinal()]);
	    printf("%16s%16d\n", "nonpiped-DADs", Gstat.panhows[DADPAN.ordinal()]);
	    printf("%16s%16d\n", "nonpiped-panels", Gstat.panhows[NOPIPE.ordinal()]);

	    /* work load distribution */
	    loadmax = loadtot = Gstat.procstat[0].fcops;
	    imax = 0;
	    for (i = 1; i < procs; ++i) {
		temp = Gstat.procstat[i].fcops;
		loadtot += temp;
		if ( temp > loadmax ) {
		    loadmax = temp;
		    imax = i;
		}
	    }
	    printf("%25s%8.2f\n", "Load balance [mean/max]", loadtot/loadmax/procs);

	    /* Delays due to pipelining. */
	    waits = 0; waittime = 0;
	    for (i = 0; i < n; i += Gstat.panstat[i].size) { /* add up all panels */
		waits += Gstat.panstat[i].pipewaits;
		waittime += Gstat.panstat[i].spintime;
	    }
	    printf("%25s%8d,\tper-panel %.1f\n", "total #delays in pipeline",
		    waits, (float)waits/panels);
	    temp = waittime / procs;
	    printf("%25s%8.2f\t[%.1f%]\n", "mean spin time per-proc",
		   temp, temp/utime[FACT.ordinal()]*100);

	    /* Delays due to scheduling. */
	    waits = 0; waittime = 0;
	    for (i = 0; i < procs; ++i) {
		waits += Gstat.procstat[i].skedwaits;
		waittime += Gstat.procstat[i].skedtime;
	    }
	    printf("%25s%8d\n", "total #delays in schedule", waits);
	    temp = waittime / procs;
	    printf("%25s%8.2f\t[%.1f%]\n", "mean sched. time per-proc",
		   temp, temp/utime[FACT.ordinal()]*100);

	    /* estimated overhead in spin-locks */
	    double TMUTEX       = 2.71e-6;
	    int FLOPS_PER_LOCK  = 407;

	    cs_numbers = n + 3*supers + panels + procs;
	    itemp = cs_numbers * FLOPS_PER_LOCK;     /* translated to flops */
	    temp = (float) (cs_numbers * TMUTEX);
	    printf("mutex-lock overhead (est.) %8.2f, #locks %d, equiv. flops %e\n",
		   temp, cs_numbers, (float) itemp);
	    printf("time in critical section   %8.2f\t[%.1f%]\n",
		   cs_time/procs, cs_time/procs/utime[FACT.ordinal()]*100);

	    printf("\n---- Parallel Profile Per Panel ----\n");
	    printf("%8s%8s%16s%8s%8s%12s%8s\n", "panel", "height",
		    "factops", "[tot%]", "msec", "spin(msec)", "Mflops");
	    thresh = (float) (0.005 * loadtot);
	    loadprint = 0;
	    itemp = 0;
	    for (i = 0; i < n; i += Gstat.panstat[i].size) {
		pan = (Gstat.panstat[i]);
		if ( pan.flopcnt > thresh ) {
		    loadprint += pan.flopcnt;
		    ++itemp;
		    if ( pan.fctime != 0 ) temp = (float) (pan.flopcnt/pan.fctime*1e-6);
		    printf("%4d%4d%8d%16e%8.1f%8.2f%12.2f%8.2f\n", i, pan.size,
			    Gstat.height[i], pan.flopcnt, pan.flopcnt/loadtot*100.,
			    pan.fctime*1e3, pan.spintime*1e3, temp);
		}
	    }
	    printf("Total panels %d,  height(T) %d, height(T)/n= %.4f\n",
		   panels, Gstat.height[n], (float)Gstat.height[n]/n);
	    printf("Printed flops %e [%.1f], printed panels %d [%.1f]\n",
		    loadprint, loadprint/loadtot*100.,
		    itemp, (float)itemp/panels);

	/*    print_flops_by_height(n, panstat, height, flops_by_height);*/

	    printf("---- End ParallelProfile().\n\n");
	    fflush(stdout);
	    return 0;
	}

	/*
	 * Print the distribution of flops by the height of etree.
	 */
	static
	void
	print_flops_by_height(int n, panstat_t panstat[],
			      int height[], float flops_by_height[])
	{
	    int i, w, ht;
	    float flops;

	    for (i = 0; i < n; i += w) {
		w = panstat[i].size;
		ht = height[i];
		flops_by_height[ht] += panstat[i].flopcnt;
	    }

	    printf("\n%8s\t%8s\n", "height", "flops");
	    ht = height[n-1]; /* root */
	    for (i = 0; i <= ht; ++i) {
		flops = flops_by_height[i];
		if ( flops != 0.0 )
		    printf("%8d\t%e\n", i, flops);
	    }
	}


	/*
	 * Print the analysis of the optimal runtime.
	 */
	static
	int
	CPprofile(final int n, cp_panel_t cp_panel[], pxgstrf_shared_t pxgstrf_shared)
	{
	    Gstat_t Gstat = pxgstrf_shared.Gstat;
	    int maxpan, i, j, treecnt;
	    float eft, maxeft; /* earliest (possible) finish time */
	    float  ops[];

	    /* Find the longest (weighted) path in the elimination forest. */
	    treecnt = 0;
	    maxeft = 0;
	    for (i = Gstat.cp_firstkid[n]; i != EMPTY; i = Gstat.cp_nextkid[i]) {
	/*	printf("Root %d, height %d\n", i, height[i]);*/
		j = (pxgstrf_shared.pan_status[i].size > 0) ?
		  i : (i + pxgstrf_shared.pan_status[i].size);
		eft   = cp_panel[j].est + cp_panel[j].pdiv;
		if ( eft > maxeft ) {
		    maxeft = eft;
		    maxpan = j;
		}
		++treecnt;
	    }

	    ops   = Gstat.ops;
	    printf("\n** Runtime prediction model: #trees %d\n", treecnt);
	    printf("Last panel %d, seq-time %e, EFT %e, ideal-speedup %.2f\n",
		   maxpan, ops[FACT.ordinal()], maxeft, ops[FACT.ordinal()]/maxeft);

	if ( DEBUGlevel>=2 ) {
	    printf("last panel %d\n", maxpan);
	    for (i = 0; i < n; i += pxgstrf_shared.pan_status[i].size)
		printf("%6d %8s%e\t%8s%8.0f\n", i, "est  ", cp_panel[i].est,
		       "pdiv  ", cp_panel[i].pdiv);
	}
	    return 0;
	}


	/***************************************************************
	 * Utilities to print the supermatrix.
	 ***************************************************************/

	static int PERLINE = 10;
	static String FMT      = "%7.4f ";

	/* A is of type Stype==SCP */
	static
	void
	Print_SuperNode_SCP(SuperMatrix A)
	{
	    int i, j, c;
	    int n = A.ncol;
	    SCPformat Astore = (SCPformat) A.Store;
	    double nzval[] = Astore.nzval;
	    int colbeg[] = Astore.nzval_colbeg, colend[] = Astore.nzval_colend;
	    printf("SuperNode_SCP: nnz %d, nsuper %d\n", Astore.nnz, Astore.nsuper);
	    printf("valL=[\n");
	    for (c = 0, j = 0; j < n; ++j) {
	        for (i = colbeg[j]; i < colend[j]; ++i) {
		    if (c == PERLINE) { printf("\n"); c = 0; }
		    printf(FMT, nzval[i]);
		    ++c;
		}
	    }
	    printf("];\n");
	    fflush(stdout);
	    /*    SUPERLU_FREE ( ((SCPformat *)A.Store).rowind );
	    SUPERLU_FREE ( ((SCPformat *)A.Store).rowind_colbeg );
	    SUPERLU_FREE ( ((SCPformat *)A.Store).rowind_colend );
	    SUPERLU_FREE ( ((SCPformat *)A.Store).col_to_sup );
	    SUPERLU_FREE ( ((SCPformat *)A.Store).sup_to_colbeg );
	    SUPERLU_FREE ( ((SCPformat *)A.Store).sup_to_colend );*/
	}

	/* A is of type NCP */
	static
	void
	Print_CompCol_NC(SuperMatrix A)
	{
	    int i, j, c;
	    int n = A.ncol;
	    NCformat Astore = (NCformat) A.Store;
	    double nzval[] = Astore.nzval;
	    int colptr[] = Astore.colptr;
	    printf("CompCol_NC: nnz %d\n", Astore.nnz);
	    printf("valA=[\n");
	    for (c = 0, j = 0; j < n; ++j) {
	        for (i = colptr[j]; i < colptr[j+1]; ++i, ++c) {
		    if (c == PERLINE) { printf("\n"); c = 0; }
		    printf(FMT, nzval[i]);
		}
	    }
	    printf("];\n");
	    fflush(stdout);
	}

	/* A is of type NCP */
	static
	void
	Print_CompCol_NCP(SuperMatrix A)
	{
	    int i, j, c;
	    int n = A.ncol;
	    NCPformat Astore = (NCPformat) A.Store;
	    double nzval[] = Astore.nzval;
	    int colbeg[] = Astore.colbeg, colend[] = Astore.colend;
	    printf("SuperNode_NCP: nnz %d\n", Astore.nnz);
	    printf("nzval[U]\n");
	    for (c = 0, j = 0; j < n; ++j) {
	        for (i = colbeg[j]; i < colend[j]; ++i, ++c) {
		    if (c == PERLINE) { printf("\n"); c = 0; }
		    printf(FMT, nzval[i]);
		}
	    }
	    printf("\n");
	    fflush(stdout);
	}

	/* A is of type DN */
	static
	void
	Print_Dense(SuperMatrix A)
	{
	    int i, j, c;
	    int m = A.nrow, n = A.ncol;
	    DNformat Astore = (DNformat) A.Store;
	    int lda = Astore.lda;
	    double nzval[] = Astore.nzval;
	    printf("Dense: lda %d\n", lda);
	    printf("val=[\n");
	    for (c = 0, j = 0; j < n; ++j) {
	        for (i = 0; i < m; ++i, ++c) {
		    if (c == PERLINE) { printf("\n"); c = 0; }
		    printf(FMT, nzval[i + j*lda]);
	      }
	    }
	    printf("];\n");
	    fflush(stdout);
	}

}

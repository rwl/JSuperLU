/*
 * -- SuperLU MT routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 */
package gov.lbl.superlu;

import gov.lbl.superlu.Dlu_pdsp_defs.GlobalLU_t;
import gov.lbl.superlu.Dlu_pxgstrf_synch.pxgstrf_relax_t;
import gov.lbl.superlu.Dlu_slu_mt_util.ExpHeader;
import gov.lbl.superlu.Dlu_slu_mt_util.MemType;
import gov.lbl.superlu.Dlu_slu_mt_util.superlu_memusage_t;
import gov.lbl.superlu.Dlu_slu_mt_util.superlumt_options_t;
import gov.lbl.superlu.Dlu_slu_mt_util.yes_no_t;
import gov.lbl.superlu.Dlu_supermatrix.NCPformat;
import gov.lbl.superlu.Dlu_supermatrix.SCPformat;
import gov.lbl.superlu.Dlu_supermatrix.SuperMatrix;

import static gov.lbl.superlu.Dlu.CHK_EXPAND;
import static gov.lbl.superlu.Dlu.DEBUG;
import static gov.lbl.superlu.Dlu.PRNTlevel;
import static gov.lbl.superlu.Dlu.fflush;
import static gov.lbl.superlu.Dlu.fprintf;
import static gov.lbl.superlu.Dlu.printf;
import static gov.lbl.superlu.Dlu.stderr;
import static gov.lbl.superlu.Dlu.stdout;
import static gov.lbl.superlu.Dlu_slu_mt_util.NO_MARKER;
import static gov.lbl.superlu.Dlu_slu_mt_util.SUPERLU_MAX;
import static gov.lbl.superlu.Dlu_slu_mt_util.MemType.LSUB;
import static gov.lbl.superlu.Dlu_slu_mt_util.MemType.LUSUP;
import static gov.lbl.superlu.Dlu_slu_mt_util.MemType.UCOL;
import static gov.lbl.superlu.Dlu_slu_mt_util.MemType.USUB;
import static gov.lbl.superlu.Dlu_slu_mt_util.yes_no_t.NO;
import static gov.lbl.superlu.Dlu_slu_mt_util.yes_no_t.YES;
import static gov.lbl.superlu.Dlu_sp_ienv.sp_ienv;



public class Dlu_pdmemory {

	/* ------------------
	   Constants & Macros
	   ------------------ */
	static final double EXPAND   = 1.5;
	static final int NO_MEMTYPE  = 4;      /* 0: lusup;
				      1: ucol;
				      2: lsub;
				      3: usub */
	static int GluIntArray(int n) {
		return 9 * (n) + 5;
	}

	static class LU_stack_t {
	    int  size;
	    int  used;
	    int  top1;  /* grow upward, relative to &array[0] */
	    int  top2;  /* grow downward */
	    Object array[];
	}

	enum stack_end_t {HEAD, TAIL}
	enum LU_space_t {SYSTEM, USER}

	static ExpHeader dexpanders[] = null; /* Array of pointers to 4 types of memory */

	static LU_stack_t stack;
	static int        no_expand;
	static int        ndim;
	static LU_space_t whichspace; /* 0 - system malloc'd; 1 - user provided */

	/* Macros to manipulate stack */
	static boolean StackFull(int x) {
		return x + stack.used >= stack.size;
	}
//	static boolean NotDoubleAlign(int addr) ( (long int)addr & 7 )
//	static boolean DoubleAlign(int addr)    ( ((long int)addr + 7) & ~7L )

	static int Reduce(int alpha) {
		return (alpha + 1) / 2;     /* i.e. (alpha-1)/2 + 1 */
	}

	/* temporary space used by BLAS calls */
	static int NUM_TEMPV(int n,int w,int t,int b) {
		return SUPERLU_MAX( 2*n, (t + b)*w );
	}

	/*
	 * Setup the memory model to be used for factorization.
	 *    lwork = 0: use system malloc;
	 *    lwork > 0: use user-supplied work[] space.
	 */
	static
	void pdgstrf_SetupSpace(Object work[], int lwork)
	{
	    if ( lwork == 0 ) {
	        whichspace = LU_space_t.SYSTEM; /* malloc/free */
	    } else if ( lwork > 0 ) {
	        whichspace = LU_space_t.USER;   /* user provided space */
	        stack.size = lwork;
	        stack.used = 0;
	        stack.top1 = 0;
	        stack.top2 = lwork;
	        stack.array = (Object []) work;
	    }
	}


	static
	Object[] duser_malloc(int bytes, int which_end)
	{
	    Object buf[];

	    if ( StackFull(bytes) ) return (null);

	    if ( which_end == HEAD ) {
		buf = (char[]) stack.array + stack.top1;
		stack.top1 += bytes;
	    } else {
		stack.top2 -= bytes;
		buf = (char[]) stack.array + stack.top2;
	    }

	    stack.used += bytes;
	    return buf;
	}


	static
	void duser_free(int bytes, int which_end)
	{
	    if ( which_end == stack_end_t.HEAD.ordinal() ) {
		stack.top1 -= bytes;
	    } else {
		stack.top2 += bytes;
	    }
	    stack.used -= bytes;
	}


	/* Returns the working storage used during factorization */
	static
	int superlu_dTempSpace(int n, int w, int p)
	{
	    float tmp, ptmp;
	    int iword = 32/*sizeof(int)*/, dword = 64/*sizeof(double)*/;
	    int    maxsuper = sp_ienv(3),
	           rowblk   = sp_ienv(4);

	    /* globally shared */
	    tmp = 14 * n * iword;

	    /* local to each processor */
	    ptmp = (2 * w + 5 + NO_MARKER) * n * iword;
	    ptmp += (n * w + NUM_TEMPV(n,w,maxsuper,rowblk)) * dword;
	if ( PRNTlevel>=1 ) {
	    printf("Per-processor work[] %.0f MB\n", ptmp/1024/1024);
	}
	    ptmp *= p;

	    return (tmp + ptmp);
	}

	/*
	 * superlu_memusage consists of the following fields:
	 *    o for_lu (float)
	 *      The amount of space used in bytes for L\U data structures.
	 *    o total_needed (float)
	 *      The amount of space needed in bytes to perform factorization.
	 *    o expansions (int)
	 *      The number of memory expansions during the LU factorization.
	 */
	static
	int superlu_dQuerySpace(int P, SuperMatrix L, SuperMatrix U, int panel_size,
	                       superlu_memusage_t superlu_memusage)
	{
	    SCPformat Lstore;
	    NCPformat Ustore;
	    int n, iword, dword, lwork;

	    Lstore = (SCPformat) L.Store;
	    Ustore = (NCPformat) U.Store;
	    n = L.ncol;
	    iword = 32/*sizeof(int)*/;
	    dword = 64/*sizeof(double)*/;

	    /* L supernodes of type SCP */
	    superlu_memusage.for_lu = (float) (7*n + 3) * iword
	                             + (float) Lstore.nzval_colend[n-1] * dword
	                             + (float) Lstore.rowind_colend[n-1] * iword;

	    /* U columns of type NCP */
	    superlu_memusage.for_lu += (2*n + 1) * iword
	        + (float) Ustore.colend[n-1] * (dword + iword);

	    /* Working storage to support factorization */
	    lwork = superlu_dTempSpace(n, panel_size, P);
	    superlu_memusage.total_needed = superlu_memusage.for_lu + lwork;

	    superlu_memusage.expansions = --no_expand;

	    return 0;
	}


	static
	float pdgstrf_memory_use(final int nzlmax, final int nzumax, final int nzlumax)
	{
	    float iword, dword, t;

	    iword   = 32/*sizeof(int)*/;
	    dword   = 64/*sizeof(double)*/;

	    t = 10.f * ndim * iword + nzlmax * iword + nzumax * (iword + dword)
		+ nzlumax * dword;
	    return t;
	}


	/*
	 * Allocate storage for the data structures common to all factor routines.
	 * For those unpredictable size, make a guess as FILL * nnz(A).
	 * Return value:
	 *     If lwork = -1, return the estimated amount of space required;
	 *     otherwise, return the amount of space actually allocated when
	 *     memory allocation failure occurred.
	 */
	static
	float
	pdgstrf_MemInit(int n, int annz, superlumt_options_t superlumt_options,
			SuperMatrix L, SuperMatrix U, GlobalLU_t Glu)
	{
	    int nprocs = superlumt_options.nprocs;
	    yes_no_t refact = superlumt_options.refact;
	    int panel_size = superlumt_options.panel_size;
	    int lwork = superlumt_options.lwork;
	    double     work[] = superlumt_options.work;
	    int      iword, dword, retries = 0;
	    SCPformat Lstore;
	    NCPformat Ustore;
	    int      xsup[], xsup_end[], supno[];
	    int      lsub[], xlsub[], xlsub_end[];
	    double   lusup[];
	    int      xlusup[], xlusup_end[];
	    double   ucol[];
	    int      usub[], xusub[], xusub_end[];
	    int      nzlmax[] = new int[1], nzumax[] = new int[1], nzlumax[] = new int[1];
	    int      FILL_LUSUP = sp_ienv(6); /* Guess the fill-in growth for LUSUP */
	    int      FILL_UCOL = sp_ienv(7); /* Guess the fill-in growth for UCOL */
	    int      FILL_LSUB = sp_ienv(8); /* Guess the fill-in growth for LSUB */

	    no_expand = 0;
	    ndim      = n;
	    iword     = 32/*sizeof(int)*/;
	    dword     = 64/*sizeof(double)*/;

	    if ( dexpanders == null )
	      dexpanders = new ExpHeader[NO_MEMTYPE];

	    if ( refact == NO ) {

		/* Guess amount of storage needed by L\U factors. */
	    if ( FILL_UCOL < 0 ) nzumax[0] = -FILL_UCOL * annz;
		else nzumax[0] = FILL_UCOL;
		if ( FILL_LSUB < 0 ) nzlmax[0] = -FILL_LSUB * annz;
		else nzlmax[0] = FILL_LSUB;

		if ( Glu.dynamic_snode_bound == YES.ordinal() ) {
		    if ( FILL_LUSUP < 0 ) nzlumax[0] = -FILL_LUSUP * annz;
		    else nzlumax[0] = FILL_LUSUP; /* estimate an upper bound */
		} else {
		    nzlumax[0] = Glu.nzlumax; /* preset as static upper bound */
		}

		if ( lwork == -1 ) {
		    return (GluIntArray(n) * iword +
			    superlu_dTempSpace(n, panel_size, nprocs)
			    + (nzlmax[0]+nzumax[0])*iword + (nzlumax[0]+nzumax[0])*dword);
	        } else {
		    pdgstrf_SetupSpace(work, lwork);
		}

		/* Integer pointers for L\U factors */
		if ( whichspace == LU_space_t.SYSTEM ) {
		    xsup       = intMalloc(n+1);
		    xsup_end   = intMalloc(n);
		    supno      = intMalloc(n+1);
		    xlsub      = intMalloc(n+1);
		    xlsub_end  = intMalloc(n);
		    xlusup     = intMalloc(n+1);
		    xlusup_end = intMalloc(n);
		    xusub      = intMalloc(n+1);
		    xusub_end  = intMalloc(n);
		} else {
		    xsup       = (int [])duser_malloc((n+1) * iword, stack_end_t.HEAD);
		    xsup_end   = (int [])duser_malloc((n) * iword, stack_end_t.HEAD);
		    supno      = (int [])duser_malloc((n+1) * iword, stack_end_t.HEAD);
		    xlsub      = (int [])duser_malloc((n+1) * iword, stack_end_t.HEAD);
		    xlsub_end  = (int [])duser_malloc((n) * iword, stack_end_t.HEAD);
		    xlusup     = (int [])duser_malloc((n+1) * iword, stack_end_t.HEAD);
		    xlusup_end = (int [])duser_malloc((n) * iword, stack_end_t.HEAD);
		    xusub      = (int [])duser_malloc((n+1) * iword, stack_end_t.HEAD);
		    xusub_end  = (int [])duser_malloc((n) * iword, stack_end_t.HEAD);
		}

		lusup = (double []) pdgstrf_expand( nzlumax, LUSUP, 0, 0, Glu );
		ucol  = (double []) pdgstrf_expand( nzumax, UCOL, 0, 0, Glu );
		lsub  = (int [])    pdgstrf_expand( nzlmax, LSUB, 0, 0, Glu );
		usub  = (int [])    pdgstrf_expand( nzumax, USUB, 0, 1, Glu );

		while ( ucol == null || lsub == null || usub == null ) {
		    /*SUPERLU_ABORT("Not enough core in LUMemInit()");*/
	if (PRNTlevel==1) {
		    printf(".. pdgstrf_MemInit(): #retries %d\n", ++retries);
	}
		    if ( whichspace == LU_space_t.SYSTEM ) {
			ucol = null;
			lsub = null;
			usub = null;
		    } else {
			duser_free(nzumax[0]*dword+(nzlmax[0]+nzumax[0])*iword, stack_end_t.HEAD);
		    }
		    nzumax[0] /= 2;    /* reduce request */
		    nzlmax[0] /= 2;
		    if ( nzumax[0] < annz/2 ) {
			printf("Not enough memory to perform factorization.\n");
			return (pdgstrf_memory_use(nzlmax[0], nzumax[0], nzlumax[0]) + n);
		    }
		    ucol  = (double []) pdgstrf_expand( nzumax, UCOL, 0, 0, Glu );
		    lsub  = (int [])    pdgstrf_expand( nzlmax, LSUB, 0, 0, Glu );
		    usub  = (int [])    pdgstrf_expand( nzumax, USUB, 0, 1, Glu );
		}

		if ( lusup == null )  {
		    float t = pdgstrf_memory_use(nzlmax[0], nzumax[0], nzlumax[0]) + n;
		    printf("Not enough memory to perform factorization .. " +
			   "need %.1f GBytes\n", t*1e-9);
		    fflush(stdout);
		    return (t);
		}

	    } else { /* refact == YES */
		Lstore   = (SCPformat) L.Store;
		Ustore   = (NCPformat) U.Store;
		xsup     = Lstore.sup_to_colbeg;
		xsup_end = Lstore.sup_to_colend;
		supno    = Lstore.col_to_sup;
		xlsub    = Lstore.rowind_colbeg;
		xlsub_end= Lstore.rowind_colend;
		xlusup   = Lstore.nzval_colbeg;
		xlusup_end= Lstore.nzval_colend;
		xusub    = Ustore.colbeg;
		xusub_end= Ustore.colend;
		nzlmax[0]   = Glu.nzlmax;    /* max from previous factorization */
		nzumax[0]   = Glu.nzumax;
		nzlumax[0]  = Glu.nzlumax;

		if ( lwork == -1 ) {
		    return (GluIntArray(n) * iword + superlu_dTempSpace(n, panel_size, nprocs)
			    + (nzlmax[0]+nzumax[0])*iword + (nzlumax[0]+nzumax[0])*dword);
	        } else if ( lwork == 0 ) {
		    whichspace = LU_space_t.SYSTEM;
		} else {
		    whichspace = LU_space_t.USER;
		    stack.size = lwork;
		    stack.top2 = lwork;
		}

		lsub  = dexpanders[LSUB.ordinal()].mem  = Lstore.rowind;
		lusup = dexpanders[LUSUP.ordinal()].mem = Lstore.nzval;
		usub  = dexpanders[USUB.ordinal()].mem  = Ustore.rowind;
		ucol  = dexpanders[UCOL.ordinal()].mem  = Ustore.nzval;;

		dexpanders[LSUB.ordinal()].size         = nzlmax[0];
		dexpanders[LUSUP.ordinal()].size        = nzlumax[0];
		dexpanders[USUB.ordinal()].size         = nzumax[0];
		dexpanders[UCOL.ordinal()].size         = nzumax[0];
	    }

	    Glu.xsup       = xsup;
	    Glu.xsup_end   = xsup_end;
	    Glu.supno      = supno;
	    Glu.lsub       = lsub;
	    Glu.xlsub      = xlsub;
	    Glu.xlsub_end  = xlsub_end;
	    Glu.lusup      = lusup;
	    Glu.xlusup     = xlusup;
	    Glu.xlusup_end = xlusup_end;
	    Glu.ucol       = ucol;
	    Glu.usub       = usub;
	    Glu.xusub      = xusub;
	    Glu.xusub_end  = xusub_end;
	    Glu.nzlmax     = nzlmax[0];
	    Glu.nzumax     = nzumax[0];
	    Glu.nzlumax    = nzlumax[0];
	    ++no_expand;

	if ( PRNTlevel>=1 ) {
	    printf(".. pdgstrf_MemInit() refact %d, space? %d, nzlumax %d, nzumax %d, nzlmax %d\n",
		refact, whichspace, nzlumax, nzumax, nzlmax);
	    printf(".. pdgstrf_MemInit() FILL_LUSUP %d, FILL_UCOL %d, FILL_LSUB %d\n",
		FILL_LUSUP, FILL_UCOL, FILL_LSUB);
	    fflush(stdout);
	}

	    return 0;

	} /* pdgstrf_MemInit */

	/*
	 * Allocate known working storage. Returns 0 if success, otherwise
	 * returns the number of bytes allocated so far when failure occurred.
	 */
	static
	int
	pdgstrf_WorkInit(int n, int panel_size, int iworkptr[][], double dworkptr[][])
	{
	    int  isize, dsize, extra;
	    double old_ptr[];
	    int    maxsuper = sp_ienv(3),
	           rowblk   = sp_ienv(4);

	    isize = (2*panel_size + 5 + NO_MARKER) * n * 32/*sizeof(int)*/;
	    dsize = (n * panel_size +
		     NUM_TEMPV(n,panel_size,maxsuper,rowblk)) * 64/*sizeof(double)*/;

	    if ( whichspace == LU_space_t.SYSTEM )
		iworkptr[0] = (int []) intCalloc(isize/32/*sizeof(int)*/);
	    else
		iworkptr[0] = (int []) duser_malloc(isize, stack_end_t.TAIL);
	    if ( iworkptr[0] == null ) {
		fprintf(stderr, "pdgstrf_WorkInit: malloc fails for local iworkptr[]\n");
		return (isize + n);
	    }

	    if ( whichspace == LU_space_t.SYSTEM )
		dworkptr[0] = (double []) SUPERLU_MALLOC(dsize);
	    else {
		dworkptr[0] = (double []) duser_malloc(dsize, stack_end_t.TAIL);
		if ( NotDoubleAlign(dworkptr[0]) ) {
		    old_ptr = dworkptr[0];
		    dworkptr[0] = (double[]) DoubleAlign(dworkptr[0]);
		    dworkptr[0] = (double[]) ((double[])dworkptr[0] - 1);
		    extra = (char[])old_ptr - (char[])dworkptr[0];
	if (CHK_EXPAND) {
		    printf("pdgstrf_WorkInit: not aligned, extra %d\n", extra);
	}
		    stack.top2 -= extra;
		    stack.used += extra;
		}
	    }
	    if ( dworkptr[0] == null ) {
		fprintf(stderr, "malloc fails for local dworkptr[].");
		return (isize + dsize + n);
	    }

	    return 0;
	}


	/*
	 * Set up pointers for real working arrays.
	 */
	static
	void
	pdgstrf_SetRWork(int n, int panel_size, double dworkptr[],
			 double dense[][], double tempv[][])
	{
	    double zero = 0.0;

	    int maxsuper = sp_ienv(3);
	    int rowblk   = sp_ienv(4);
	    dense[0] = dworkptr;
	    tempv[0] = dense[0] + panel_size*n;
	    dfill (dense[0], n * panel_size, zero);
	    dfill (tempv[0], NUM_TEMPV(n,panel_size,maxsuper,rowblk), zero);
	}

	/*
	 * Free the working storage used by factor routines.
	 */
	static
	void pdgstrf_WorkFree(int iwork[], double dwork[], GlobalLU_t Glu)
	{
	    if ( whichspace == LU_space_t.SYSTEM ) {
		iwork = null;
		dwork = null;
	    } else {
		stack.used -= (stack.size - stack.top2);
		stack.top2 = stack.size;
	/*	pdgstrf_StackCompress(Glu);  */
	    }
	}

	/*
	 * Expand the data structures for L and U during the factorization.
	 * Return value:   0 - successful return
	 *               > 0 - number of bytes allocated when run out of space
	 */
	int
	pdgstrf_MemXpand(
			 int jcol,
			 int next, /* number of elements currently in the factors */
			 MemType mem_type,/* which type of memory to expand  */
			 int maxlen[], /* modified - max. length of a data structure */
			 GlobalLU_t Glu /* modified - global LU data structures */
			 )
	{
	    Object   new_mem[];

	if (CHK_EXPAND) {
	    printf("pdgstrf_MemXpand(): jcol %d, next %d, maxlen %d, MemType %d\n",
		   jcol, next, maxlen[0], mem_type);
	}

	    if (mem_type == USUB)
	    	new_mem = pdgstrf_expand(maxlen, mem_type, next, 1, Glu);
	    else
		new_mem = pdgstrf_expand(maxlen, mem_type, next, 0, Glu);

	    if ( new_mem == null ) {
		int    nzlmax  = Glu.nzlmax;
		int    nzumax  = Glu.nzumax;
		int    nzlumax = Glu.nzlumax;
	    	fprintf(stderr, "Can't expand MemType %d: jcol %d\n", mem_type, jcol);
	    	return (pdgstrf_memory_use(nzlmax, nzumax, nzlumax) + ndim);
	    }

	    switch ( mem_type ) {
	      case LUSUP:
		Glu.lusup   = (double []) new_mem;
		Glu.nzlumax = maxlen;
		break;
	      case UCOL:
		Glu.ucol   = (double []) new_mem;
		Glu.nzumax = maxlen;
		break;
	      case LSUB:
		Glu.lsub   = (int []) new_mem;
		Glu.nzlmax = maxlen;
		break;
	      case USUB:
		Glu.usub   = (int []) new_mem;
		Glu.nzumax = maxlen;
		break;
	    }

	    return 0;

	}


	static
	void
	copy_mem_double(int howmany, double old[], double new_[])
	{
	    int i;
	    double dold[] = old;
	    double dnew[] = new_;
	    for (i = 0; i < howmany; i++) dnew[i] = dold[i];
	}


	/*
	 * Expand the existing storage to accommodate more fill-ins.
	 */
	static
	Object[]
	pdgstrf_expand(
	                int prev_len[],  /* length used from previous call */
	                MemType type,    /* which part of the memory to expand */
	                int len_to_copy, /* size of memory to be copied to new store */
	                int keep_prev,   /* = 1: use prev_len;
	                                    = 0: compute new_len to expand */
	                GlobalLU_t Glu   /* modified - global LU data structures */
	                )
	{
	    double   alpha = EXPAND;
	    Object   new_mem[], old_mem[];
	    int      new_len, tries, lword, extra, bytes_to_copy;

	    if ( no_expand == 0 || keep_prev != 0 ) /* First time allocate requested */
	        new_len = prev_len[0];
	    else {
	        new_len = (int) (alpha * prev_len[0]);
	    }

	    if ( type == LSUB || type == USUB ) lword = 32/*sizeof(int)*/;
	    else lword = 64/*sizeof(double)*/;

	    if ( whichspace == LU_space_t.SYSTEM ) {
	        new_mem = (Object []) SUPERLU_MALLOC( (size_t) new_len * lword );

	        if ( no_expand != 0 ) {
	            tries = 0;
	            if ( keep_prev != 0 ) {
	                if ( new_mem == null ) return (null);
	            } else {
	                while ( new_mem == null ) {
	                    if ( ++tries > 10 ) return (null);
	                    alpha = Reduce(alpha);
	                    new_len = alpha * prev_len[0];
	                    new_mem = (Object []) SUPERLU_MALLOC((size_t) new_len * lword);
	                }
	            }
	            if ( type == LSUB || type == USUB ) {
	                copy_mem_int(len_to_copy, dexpanders[type.ordinal()].mem, new_mem);
	            } else {
	                copy_mem_double(len_to_copy, dexpanders[type.ordinal()].mem, new_mem);
	            }
	            dexpanders[type.ordinal()].mem = null;
	        }
	        dexpanders[type.ordinal()].mem = (Object []) new_mem;

	    } else { /* whichspace == USER */
	        if ( no_expand == 0 ) {
	            new_mem = duser_malloc(new_len * lword, HEAD);
	            if ( NotDoubleAlign(new_mem) &&
	                (type == LUSUP || type == UCOL) ) {
	                old_mem = new_mem;
	                new_mem = (Object [])DoubleAlign(new_mem);
	                extra = (char[])new_mem - (char[])old_mem;
	if (CHK_EXPAND) {
	                printf("expand(): not aligned, extra %d\n", extra);
	}
	                stack.top1 += extra;
	                stack.used += extra;
	            }
	            dexpanders[type.ordinal()].mem = (Object []) new_mem;
	        }
	        else {
	            tries = 0;
	            extra = (new_len - prev_len[0]) * lword;
	            if ( keep_prev != 0 ) {
	                if ( StackFull(extra) ) return (null);
	            } else {
	                while ( StackFull(extra) ) {
	                    if ( ++tries > 10 ) return (null);
	                    alpha = Reduce(alpha);
	                    new_len = alpha * prev_len[0];
	                    extra = (new_len - prev_len[0]) * lword;
	                }
	            }

	            if ( type != USUB ) {
	                new_mem = (Object[])((char[])dexpanders[type.ordinal() + 1].mem + extra);
	                bytes_to_copy = (char[])stack.array + stack.top1
	                    - (char[])dexpanders[type.ordinal() + 1].mem;
	                user_bcopy(dexpanders[type.ordinal()+1].mem, new_mem, bytes_to_copy);

	                if ( type.ordinal() < USUB.ordinal() ) {
	                    Glu.usub = dexpanders[USUB.ordinal()].mem =
	                        (Object[])((char[])dexpanders[USUB.ordinal()].mem + extra);
	                }
	                if ( type.ordinal() < LSUB.ordinal() ) {
	                    Glu.lsub = dexpanders[LSUB.ordinal()].mem =
	                        (Object[])((char[])dexpanders[LSUB.ordinal()].mem + extra);
	                }
	                if ( type.ordinal() < UCOL.ordinal() ) {
	                    Glu.ucol = dexpanders[UCOL.ordinal()].mem =
	                        (Object[])((char[])dexpanders[UCOL.ordinal()].mem + extra);
	                }
	                stack.top1 += extra;
	                stack.used += extra;
	                if ( type == UCOL ) {
	                    stack.top1 += extra;   /* Add same amount for USUB */
	                    stack.used += extra;
	                }

	            } /* if ... */

	        } /* else ... */
	    }
	if (DEBUG) {
	    printf("pdgstrf_expand[type %d]\n", type.ordinal());
	}
	    dexpanders[type.ordinal()].size = new_len;
	    prev_len[0] = new_len;
	    if ( no_expand != 0 ) ++no_expand;

	    return (Object []) dexpanders[type.ordinal()].mem;

	} /* expand */


	/*
	 * Compress the work[] array to remove fragmentation.
	 */
	static
	void
	pdgstrf_StackCompress(GlobalLU_t Glu)
	{
	    int iword, dword;
	    char     last[], fragment[];
	    int      ifrom[], ito[];
	    double   dfrom[], dto[];
	    int      xlsub[], lsub[], xusub_end[], usub[], xlusup[];
	    double   ucol[], lusup[];

	    iword = 32/*sizeof(int)*/;
	    dword = 64/*sizeof(double)*/;

	    xlsub  = Glu.xlsub;
	    lsub   = Glu.lsub;
	    xusub_end  = Glu.xusub_end;
	    usub   = Glu.usub;
	    xlusup = Glu.xlusup;
	    ucol   = Glu.ucol;
	    lusup  = Glu.lusup;

	    dfrom = ucol;
	    dto = (double [])((char[])lusup + xlusup[ndim] * dword);
	    copy_mem_double(xusub_end[ndim-1], dfrom, dto);
	    ucol = dto;

	    ifrom = lsub;
	    ito = (int []) ((char[])ucol + xusub_end[ndim-1] * iword);
	    copy_mem_int(xlsub[ndim], ifrom, ito);
	    lsub = ito;

	    ifrom = usub;
	    ito = (int []) ((char[])lsub + xlsub[ndim] * iword);
	    copy_mem_int(xusub_end[ndim-1], ifrom, ito);
	    usub = ito;

	    last = (char[])usub + xusub_end[ndim-1] * iword;
	    fragment = (char[]) ((char[])stack.array + stack.top1 - last);
	    stack.used -= (long) fragment;
	    stack.top1 -= (long) fragment;

	    Glu.ucol = ucol;
	    Glu.lsub = lsub;
	    Glu.usub = usub;

	if (CHK_EXPAND) {
	    printf("pdgstrf_StackCompress: fragment %d\n", fragment);
	    /* PrintStack("After compress", Glu);
	    for (last = 0; last < ndim; ++last)
		print_lu_col("After compress:", last, 0);*/
	}

	}

	/*
	 * Allocate storage for original matrix A
	 */
	static
	void
	dallocateA(int n, int nnz, double a[][], int asub[][], int xa[][])
	{
	    a[0]    = (double []) doubleMalloc(nnz);
	    asub[0] = (int []) intMalloc(nnz);
	    xa[0]   = (int []) intMalloc(n+1);
	}

	static
	double[] doubleMalloc(int n)
	{
	    double buf[];
	    buf = (double []) SUPERLU_MALLOC( (size_t) n * 64/*sizeof(double)*/ );
	    if ( buf == null ) {
		fprintf(stderr, "SUPERLU_MALLOC failed for buf in doubleMalloc()");
		exit (1);
	    }
	    return (buf);
	}

	static
	double[] doubleCalloc(int n)
	{
	    double[] buf;
	    int i;
	    double zero = 0.0;
	    buf = (double []) SUPERLU_MALLOC( (size_t) n * 64/*sizeof(double)*/ );
	    if ( !buf ) {
		fprintf(stderr, "SUPERLU_MALLOC failed for buf in doubleCalloc()");
		exit (1);
	    }
	    for (i = 0; i < n; ++i) buf[i] = zero;
	    return (buf);
	}

	/*
	 * Set up memory image in lusup[*], using the supernode boundaries in
	 * the Householder matrix.
	 *
	 * In both static and dynamic scheme, the relaxed supernodes (leaves)
	 * are stored in the beginning of lusup[*]. In the static scheme, the
	 * memory is also set aside for the internal supernodes using upper
	 * bound information from H. In the dynamic scheme, however, the memory
	 * for the internal supernodes is not allocated by this routine.
	 *
	 * Return value
	 *   o Static scheme: number of nonzeros of all the supernodes in H.
	 *   o Dynamic scheme: number of nonzeros of the relaxed supernodes.
	 */
	static
	int
	dPresetMap(
		  final int n,
		  SuperMatrix A, /* original matrix permuted by columns */
		  pxgstrf_relax_t pxgstrf_relax[], /* relaxed supernodes */
		  superlumt_options_t superlumt_options, /* input */
		  GlobalLU_t Glu /* modified */
		  )
	{
	    int i, j, k, w, rs, rs_lastcol, krow, kmark, maxsup, nextpos;
	    int rs_nrow; /* number of nonzero rows in a relaxed supernode */
	    int          marker[], asub[], xa_begin[], xa_end[];
	    NCPformat    Astore;
	    int map_in_sup[]; /* memory mapping function; values irrelevant on entry. */
	    int colcnt[];     /* column count of Lc or H */
	    int super_bnd[];  /* supernodes partition in H */
	    char snode_env[]/*, *getenv()*/;

	    snode_env = getenv("SuperLU_DYNAMIC_SNODE_STORE");
	    if ( snode_env != null ) {
		Glu.dynamic_snode_bound = YES;
	if ( PRNTlevel>=1 ) {
		printf(".. Use dynamic alg. to allocate storage for L supernodes.\n");
	}
	    } else  Glu.dynamic_snode_bound = NO;

	    Astore   = A.Store;
	    asub     = Astore.rowind;
	    xa_begin = Astore.colbeg;
	    xa_end   = Astore.colend;
	    rs       = 1;
	    marker   = intMalloc(n);
	    ifill(marker, n, EMPTY);
	    map_in_sup = Glu.map_in_sup = intCalloc(n+1);
	    colcnt = superlumt_options.colcnt_h;
	    super_bnd = superlumt_options.part_super_h;
	    nextpos = 0;

	    /* Split large supernode into smaller pieces */
	    maxsup = sp_ienv(3);
	    for (j = 0; j < n; ) {
		w = super_bnd[j];
		k = j + w;
		if ( w > maxsup ) {
		    w = w % maxsup;
		    if ( w == 0 ) w = maxsup;
		    while ( j < k ) {
			super_bnd[j] = w;
			j += w;
			w = maxsup;
		    }
		}
		j = k;
	    }

	    for (j = 0; j < n; j += w) {
	        if ( Glu.dynamic_snode_bound == NO ) map_in_sup[j] = nextpos;

		if ( pxgstrf_relax[rs].fcol == j ) {
		    /* Column j starts a relaxed supernode. */
		    map_in_sup[j] = nextpos;
		    rs_nrow = 0;
		    w = pxgstrf_relax[rs++].size;
		    rs_lastcol = j + w;
		    for (i = j; i < rs_lastcol; ++i) {
			/* for each nonzero in A[*,i] */
			for (k = xa_begin[i]; k < xa_end[i]; k++) {
			    krow = asub[k];
			    kmark = marker[krow];
			    if ( kmark != j ) { /* first time visit krow */
				marker[krow] = j;
				++rs_nrow;
			    }
			}
		    }
		    nextpos += w * rs_nrow;

		    /* Find the next H-supernode, with leading column i, which is
		       outside the relaxed supernode, rs. */
		    for (i = j; i < rs_lastcol; k = i, i += super_bnd[i]);
		    if ( i > rs_lastcol ) {
			/* The w columns [rs_lastcol, i) may join in the
			   preceeding relaxed supernode; make sure we leave
			   enough room for the combined supernode. */
			w = i - rs_lastcol;
			nextpos += w * SUPERLU_MAX( rs_nrow, colcnt[k] );
		    }
		    w = i - j;
		} else { /* Column j starts a supernode in H */
		    w = super_bnd[j];
		    if ( Glu.dynamic_snode_bound == NO ) nextpos += w * colcnt[j];
		}

		/* Set up the offset (negative) to the leading column j of a
		   supernode in H */
		for (i = 1; i < w; ++i) map_in_sup[j + i] = -i;

	    } /* for j ... */

	    if ( Glu.dynamic_snode_bound == YES ) Glu.nextlu = nextpos;
	    else map_in_sup[n] = nextpos;

	if ( PRNTlevel>=1 ) {
	    printf("** PresetMap() allocates %d reals to lusup[*]....\n", nextpos);
	}

	    free (marker);
	    return nextpos;
	}

}

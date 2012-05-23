/*! @file Dlu_dmemory.java
 * \brief Memory details
 *
 * <pre>
 * -- SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Laboratory.
 * June 30, 2009
 * </pre>
 */
package lu.jsuper;

import lu.jsuper.Dlu_slu_ddefs.GlobalLU_t;
import lu.jsuper.Dlu_superlu_enum_consts.LU_space_t;
import lu.jsuper.Dlu_superlu_enum_consts.MemType;
import lu.jsuper.Dlu_superlu_enum_consts.fact_t;
import lu.jsuper.Dlu_superlu_enum_consts.stack_end_t;
import lu.jsuper.Dlu_supermatrix.NCformat;
import lu.jsuper.Dlu_supermatrix.SCformat;
import lu.jsuper.Dlu_supermatrix.SuperMatrix;

import static lu.jsuper.Dlu_slu_util.ABORT;
import static lu.jsuper.Dlu_slu_util.ExpHeader;
import static lu.jsuper.Dlu_slu_util.NO_MEMTYPE;
import static lu.jsuper.Dlu_slu_util.SUPERLU_MAX;
import static lu.jsuper.Dlu_slu_util.GluIntArray;
import static lu.jsuper.Dlu_slu_util.NO_MARKER;
import static lu.jsuper.Dlu_slu_util.NUM_TEMPV;

import static lu.jsuper.Dlu_util.PRNTlevel;
import static lu.jsuper.Dlu_util.DEBUG;

import static lu.jsuper.Dlu_memory.intMalloc;
import static lu.jsuper.Dlu_memory.copy_mem_int;

import static lu.jsuper.Dlu_sp_ienv.sp_ienv;

import static lu.jsuper.Dlu_dutil.dfill;


public class Dlu_dmemory {

	public static int StackFull(int x, GlobalLU_t Glu) {
		return x + (Glu.stack.used >= Glu.stack.size ? 1 : 0);
	}

	public static int NotDoubleAlign(int addr) {
		return addr & 7;
	}

	public static int TempSpace(int m, int w) {
		return (2*w + 4 + NO_MARKER) * m * 32 +//sizeof(int) +
		      (w + 1) * m * 64;//sizeof(double)
	}

	public static float Reduce(float alpha) {
		return (alpha + 1) / 2;  /* i.e. (alpha-1)/2 + 1 */
	}

	/*! \brief Setup the memory model to be used for factorization.
	 *
	 *    lwork = 0: use system malloc;
	 *    lwork > 0: use user-supplied work[] space.
	 */
	public static void dSetupSpace(double work[], int lwork, GlobalLU_t Glu) {
	    if ( lwork == 0 ) {
		Glu.MemModel = LU_space_t.SYSTEM; /* malloc/free */
	    } else if ( lwork > 0 ) {
		Glu.MemModel = LU_space_t.USER;   /* user provided space */
		Glu.stack.used = 0;
		Glu.stack.top1 = 0;
		Glu.stack.top2 = (lwork/4)*4; /* must be word addressable */
		Glu.stack.size = Glu.stack.top2;
		Glu.stack.array = (double[]) work;
	    }
	}

	public static double[] duser_malloc(int bytes, int which_end, GlobalLU_t Glu)
	{
	    double[] buf;
	    int buf_offset;

	    if ( StackFull(bytes, Glu) != 0 ) return (null);

	    if ( which_end == stack_end_t.HEAD.ordinal() ) {
		buf = Glu.stack.array;
		buf_offset = Glu.stack.top1;
		Glu.stack.top1 += bytes;
	    } else {
		Glu.stack.top2 -= bytes;
		buf = Glu.stack.array;
		buf_offset = Glu.stack.top2;
	    }

	    Glu.stack.used += bytes;
	    return buf;
	}

	/**! \brief Allocate storage for the data structures common to all factor routines.
	 *
	 * <pre>
	 * For those unpredictable size, estimate as fill_ratio * nnz(A).
	 * Return value:
	 *     If lwork = -1, return the esti()mated amount of space required, plus n;
	 *     otherwise, return the amount of space actually allocated when
	 *     memory allocation failure occurred.
	 * </pre>
	 */
	public static int dLUMemInit(fact_t fact, double work[], int lwork,
			int m, int n, int annz, int panel_size, double fill_ratio,
			SuperMatrix L, SuperMatrix U, GlobalLU_t Glu, int iwork[][],
			double dwork[][]) {

	    int      info, iword, dword;
	    SCformat Lstore;
	    NCformat Ustore;
	    int      xsup[], supno[];
	    int      lsub[], xlsub[];
	    double   lusup[];
	    int      xlusup[];
	    double   ucol[];
	    int      usub[], xusub[];
	    int      nzlmax, nzumax, nzlumax;

	    iword     = 32; //sizeof(int);
	    dword     = 64; //sizeof(double);
	    Glu.n     = n;
	    Glu.num_expansions = 0;

	    if ( Glu.expanders == null )
	        Glu.expanders = new ExpHeader[NO_MEMTYPE];

	    if ( fact != fact_t.SamePattern_SameRowPerm ) {
    	/* Guess for L\U factors */
    	nzumax = nzlumax = (int) fill_ratio * annz;
    	nzlmax = (int) (SUPERLU_MAX(1, fill_ratio/4.) * annz);

    	if ( lwork == -1 ) {
    	    return ( GluIntArray(n) * iword + TempSpace(m, panel_size)
    		    + (nzlmax+nzumax)*iword + (nzlumax+nzumax)*dword + n );
            } else {
    	    dSetupSpace(work, lwork, Glu);
    	}

    	if ( PRNTlevel >= 1 ) {
    		System.out.printf("dLUMemInit() called: fill_ratio %.0f, nzlmax %ld, nzumax %ld\n",
    		       fill_ratio, nzlmax, nzumax);
    		System.out.flush();
    	}

    	/* Integer pointers for L\U factors */
    	if ( Glu.MemModel == LU_space_t.SYSTEM ) {
    	    xsup   = intMalloc(n+1);
    	    supno  = intMalloc(n+1);
    	    xlsub  = intMalloc(n+1);
    	    xlusup = intMalloc(n+1);
    	    xusub  = intMalloc(n+1);
    	} else {
    		throw new UnsupportedOperationException();
//    	    xsup   = (int [])duser_malloc((n+1) * iword, stack_end_t.HEAD.ordinal(), Glu);
//    	    supno  = (int [])duser_malloc((n+1) * iword, stack_end_t.HEAD.ordinal(), Glu);
//    	    xlsub  = (int [])duser_malloc((n+1) * iword, stack_end_t.HEAD.ordinal(), Glu);
//    	    xlusup = (int [])duser_malloc((n+1) * iword, stack_end_t.HEAD.ordinal(), Glu);
//    	    xusub  = (int [])duser_malloc((n+1) * iword, stack_end_t.HEAD.ordinal(), Glu);
    	}

    	lusup = (double *) dexpand( nzlumax, MemType.LUSUP, 0, 0, Glu );
    	ucol  = (double *) dexpand( nzumax, MemType.UCOL, 0, 0, Glu );
    	lsub  = (int *)    dexpand( nzlmax, MemType.LSUB, 0, 0, Glu );
    	usub  = (int *)    dexpand( nzumax, MemType.USUB, 0, 1, Glu );

    	while ( lusup == 0 || ucol == 0 || lsub == 0 || usub == 0 ) {
    	    if ( Glu.MemModel == LU_space_t.SYSTEM ) {
    		SUPERLU_FREE(lusup);
    		SUPERLU_FREE(ucol);
    		SUPERLU_FREE(lsub);
    		SUPERLU_FREE(usub);
    	    } else {
    		duser_free((nzlumax+nzumax)*dword+(nzlmax+nzumax)*iword,
                                stack_end_t.HEAD, Glu);
    	    }
    	    nzlumax /= 2;
    	    nzumax /= 2;
    	    nzlmax /= 2;
    	    if ( nzlumax < annz ) {
    		System.out.printf("Not enough memory to perform factorization.\n");
    		return (dmemory_usage(nzlmax, nzumax, nzlumax, n) + n);
    	    }
    	    if ( PRNTlevel >= 1) {
    	    System.out.printf("dLUMemInit() reduce size: nzlmax %ld, nzumax %ld\n",
    		   nzlmax, nzumax);
    	    System.out.flush();
    	    }
    	    lusup = (double *) dexpand( &nzlumax, LUSUP, 0, 0, Glu );
    	    ucol  = (double *) dexpand( &nzumax, UCOL, 0, 0, Glu );
    	    lsub  = (int *)    dexpand( &nzlmax, LSUB, 0, 0, Glu );
    	    usub  = (int *)    dexpand( &nzumax, USUB, 0, 1, Glu );
    	}

	    } else {
		/* fact == SamePattern_SameRowPerm */
		Lstore   = (SCformat) L.Store;
		Ustore   = (NCformat) U.Store;
		xsup     = Lstore.sup_to_col;
		supno    = Lstore.col_to_sup;
		xlsub    = Lstore.rowind_colptr;
		xlusup   = Lstore.nzval_colptr;
		xusub    = Ustore.colptr;
		nzlmax   = Glu.nzlmax;    /* max from previous factorization */
		nzumax   = Glu.nzumax;
		nzlumax  = Glu.nzlumax;

		if ( lwork == -1 ) {
		    return ( GluIntArray(n) * iword + TempSpace(m, panel_size)
			    + (nzlmax+nzumax)*iword + (nzlumax+nzumax)*dword + n );
	        } else if ( lwork == 0 ) {
		    Glu.MemModel = LU_space_t.SYSTEM;
		} else {
		    Glu.MemModel = LU_space_t.USER;
		    Glu.stack.top2 = (lwork/4)*4; /* must be word-addressable */
		    Glu.stack.size = Glu.stack.top2;
		}

		lsub  = Glu.expanders[MemType.LSUB.ordinal()].mem  = Lstore.rowind;
		lusup = Glu.expanders[MemType.LUSUP.ordinal()].mem = Lstore.nzval;
		usub  = Glu.expanders[MemType.USUB.ordinal()].mem  = Ustore.rowind;
		ucol  = Glu.expanders[MemType.UCOL.ordinal()].mem  = Ustore.nzval;;
		Glu.expanders[MemType.LSUB.ordinal()].size         = nzlmax;
		Glu.expanders[MemType.LUSUP.ordinal()].size        = nzlumax;
		Glu.expanders[MemType.USUB.ordinal()].size         = nzumax;
		Glu.expanders[MemType.UCOL.ordinal()].size         = nzumax;
	    }

	    Glu.xsup    = xsup;
	    Glu.supno   = supno;
	    Glu.lsub    = lsub;
	    Glu.xlsub   = xlsub;
	    Glu.lusup   = lusup;
	    Glu.xlusup  = xlusup;
	    Glu.ucol    = ucol;
	    Glu.usub    = usub;
	    Glu.xusub   = xusub;
	    Glu.nzlmax  = nzlmax;
	    Glu.nzumax  = nzumax;
	    Glu.nzlumax = nzlumax;

	    info = dLUWorkInit(m, n, panel_size, iwork, dwork, Glu);
	    if ( info != 0 )
		return ( info + dmemory_usage(nzlmax, nzumax, nzlumax, n) + n);

	    ++Glu.num_expansions;

	    return 0;
	}

	/*! \brief Allocate known working storage. Returns 0 if success, otherwise
	   returns the number of bytes allocated so far when failure occurred. */
	public static int dLUWorkInit(int m, int n, int panel_size, int iworkptr[][],
            double dworkptr[][], GlobalLU_t Glu) {
	    int    isize, dsize, extra;
	    double old_ptr[];
	    int    maxsuper = SUPERLU_MAX( sp_ienv(3), sp_ienv(7) ),
	           rowblk   = sp_ienv(4);

	    isize = ( (2 * panel_size + 3 + NO_MARKER ) * m + n ) * 32;//sizeof(int);
	    dsize = (m * panel_size +
		     NUM_TEMPV(m,panel_size,maxsuper,rowblk)) * 64;//sizeof(double);

	    if ( Glu.MemModel == LU_space_t.SYSTEM )
		iworkptr[0] = (int []) intCalloc(isize/32);//sizeof(int));
	    else
		iworkptr[0] = (int []) duser_malloc(isize, stack_end_t.TAIL.ordinal(), Glu);
	    if ( iworkptr[0] == null ) {
		System.err.printf("dLUWorkInit: malloc fails for local iworkptr[]\n");
		return (isize + n);
	    }

	    if ( Glu.MemModel == SYSTEM )
		dworkptr[0] = (double []) SUPERLU_MALLOC(dsize);
	    else {
		dworkptr[0] = (double []) duser_malloc(dsize, stack_end_t.TAIL.ordinal(), Glu);
		if ( NotDoubleAlign(dworkptr[0]) ) {
		    old_ptr = dworkptr[0];
		    dworkptr[0] = (double[]) DoubleAlign(dworkptr[0]);
		    dworkptr[0] = (double[]) ((double[])dworkptr[0] - 1);
//		    extra = (char*)old_ptr - (char*)*dworkptr;
		    extra = old_ptr.length - dworkptr[0].length;
		    if (DEBUG) {
		    System.out.printf("dLUWorkInit: not aligned, extra %d\n", extra);
		    }
		    Glu.stack.top2 -= extra;
		    Glu.stack.used += extra;
		}
	    }
	    if ( dworkptr[0] == null ) {
		System.err.printf("malloc fails for local dworkptr[].");
		return (isize + dsize + n);
	    }

	    return 0;
	}

	/**! \brief Set up pointers for real working arrays.
	 */
	public static void dSetRWork(int m, int panel_size, double dworkptr[],
		 double dense[][], double tempv[][]) {
	    double zero = 0.0;

	    int maxsuper = SUPERLU_MAX( sp_ienv(3), sp_ienv(7) ),
	        rowblk   = sp_ienv(4);
	    dense[0] = dworkptr;
	    tempv[0] = dense[0] + panel_size*m;
	    dfill (dense[0], m * panel_size, zero);
	    dfill (tempv[0], NUM_TEMPV(m,panel_size,maxsuper,rowblk), zero);
	}

	public static void copy_mem_double(int howmany, double old[], double new_[])
	{
	    int i;
	    double dold[] = old;
	    double dnew[] = new_;
	    for (i = 0; i < howmany; i++) dnew[i] = dold[i];
	}

    /**! \brief Expand the existing storage to accommodate more fill-ins.
     */
	public static double[] dexpand(
			 int prev_len[],  /* length used from previous call */
			 MemType type,    /* which part of the memory to expand */
			 int len_to_copy, /* size of the memory to be copied to new store */
			 int keep_prev,   /* = 1: use prev_len;
					     = 0: compute new_len to expand */
			 GlobalLU_t Glu   /* modified - global LU data structures */) {

	    float    EXPAND = 1.5f;
	    float    alpha;
	    double   new_mem[], old_mem[];
	    int      new_len, tries, lword, extra, bytes_to_copy;
	    ExpHeader expanders[] = Glu.expanders; /* Array of 4 types of memory */

	    alpha = EXPAND;

	    if ( Glu.num_expansions == 0 || keep_prev != 0 ) {
	        /* First time allocate requested */
	        new_len = prev_len[0];
	    } else {
		new_len = (int) (alpha * prev_len[0]);
	    }

	    if ( type == MemType.LSUB || type == MemType.USUB ) lword = 32;//sizeof(int);
	    else lword = 64;//sizeof(double);

	    if ( Glu.MemModel == LU_space_t.SYSTEM ) {
		new_mem = new double[new_len];// * lword];
		if ( Glu.num_expansions != 0 ) {
		    tries = 0;
		    if ( keep_prev != 0 ) {
			if ( new_mem == null ) return (null);
		    } else {
			while ( new_mem == null ) {
			    if ( ++tries > 10 ) return (null);
			    alpha = Reduce(alpha);
			    new_len = (int) alpha * prev_len[0];
			    new_mem = new double[new_len];// * lword];
			}
		    }
		    if ( type == MemType.LSUB || type == MemType.USUB ) {
			copy_mem_int(len_to_copy, expanders[type.ordinal()].mem, new_mem);
		    } else {
			copy_mem_double(len_to_copy, expanders[type.ordinal()].mem, new_mem);
		    }
		    expanders[type.ordinal()].mem = null;
		}
		expanders[type.ordinal()].mem = new_mem;

	    } else { /* MemModel == USER */
		if ( Glu.num_expansions == 0 ) {
		    new_mem = duser_malloc(new_len * lword, stack_end_t.HEAD.ordinal(), Glu);
		    if ( NotDoubleAlign(new_mem) &&
			(type == MemType.LUSUP || type == MemType.UCOL) ) {
			old_mem = new_mem;
			new_mem = (double []) DoubleAlign(new_mem);
//			extra = (char*)new_mem - (char*)old_mem;
			extra = new_mem.length - old_mem.length;

			if (DEBUG) {
			System.out.printf("expand(): not aligned, extra %d\n", extra);
			}
			Glu.stack.top1 += extra;
			Glu.stack.used += extra;
		    }
		    expanders[type.ordinal()].mem = (double []) new_mem;
		} else {
		    tries = 0;
		    extra = (new_len - prev_len[0]) * lword;
		    if ( keep_prev != 0 ) {
			if ( StackFull(extra, Glu) != 0 ) return (null);
		    } else {
			while ( StackFull(extra, Glu) != 0 ) {
			    if ( ++tries > 10 ) return (null);
			    alpha = Reduce(alpha);
			    new_len = (int) (alpha * prev_len[0]);
			    extra = (new_len - prev_len[0]) * lword;
			}
		    }

		    if ( type != MemType.USUB ) {
			new_mem = ((double[])expanders[type.ordinal() + 1].mem + extra);
			bytes_to_copy = (double[])Glu.stack.array + Glu.stack.top1
			    - (double[])expanders[type + 1].mem;
			user_bcopy(expanders[type.ordinal()+1].mem, new_mem, bytes_to_copy);

			if ( type.ordinal() < MemType.USUB.ordinal() ) {
			    Glu.usub = expanders[MemType.USUB.ordinal()].mem =
				((double[])expanders[MemType.USUB.ordinal()].mem + extra);
			}
			if ( type.ordinal() < MemType.LSUB.ordinal() ) {
			    Glu.lsub = expanders[MemType.LSUB.ordinal()].mem =
				((double[])expanders[MemType.LSUB.ordinal()].mem + extra);
			}
			if ( type.ordinal() < MemType.UCOL.ordinal() ) {
			    Glu.ucol = expanders[MemType.UCOL.ordinal()].mem =
				((double[])expanders[MemType.UCOL.ordinal()].mem + extra);
			}
			Glu.stack.top1 += extra;
			Glu.stack.used += extra;
			if ( type.ordinal() == MemType.UCOL.ordinal() ) {
			    Glu.stack.top1 += extra;   /* Add same amount for USUB */
			    Glu.stack.used += extra;
			}

		    } /* if ... */

		} /* else ... */
	    }

	    expanders[type.ordinal()].size = new_len;
	    prev_len[0] = new_len;
	    if ( Glu.num_expansions != 0 ) ++Glu.num_expansions;

	    return (double []) expanders[type.ordinal()].mem;
	}

	public static double[] doubleMalloc(int n) {
		double[] buf = null;
		try {
			buf = new double[n];
		} catch (OutOfMemoryError e) {
			ABORT("SUPERLU_MALLOC failed for buf in doubleMalloc()\n", e);
		}
		return (buf);
	}

	public static int dmemory_usage(final int nzlmax, final int nzumax,
			  final int nzlumax, final int n) {
	    int iword, dword;

	    iword   = 32;//sizeof(int);
	    dword   = 64;//sizeof(double);

	    return (10 * n * iword +
		    nzlmax * iword + nzumax * (iword + dword) + nzlumax * dword);
	}

}

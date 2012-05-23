/**! @file sp_coletree.c
 * \brief Tree layout and computation routines
 *
 *<pre>
 * -- SuperLU routine (version 3.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 1, 2008
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

/*  Elimination tree computation and layout routines */

import static lu.jsuper.Dlu_slu_util.ABORT;


public class Dlu_sp_coletree {

	public static int[] mxCallocInt(int n)
	{
	    int i;
	    int buf[] = null;

	    try {
	    	buf = new int [n];
	    } catch (OutOfMemoryError e) {
	         ABORT("SUPERLU_MALLOC fails for buf in mxCallocInt()");
	    }
	    for (i = 0; i < n; i++) buf[i] = 0;
	    return (buf);
	}

	/*
	 *  q = TreePostorder (n, p);
	 *
	 *	Postorder a tree.
	 *	Input:
	 *	  p is a vector of parent pointers for a forest whose
	 *        vertices are the integers 0 to n-1; p[root]==n.
	 *	Output:
	 *	  q is a vector indexed by 0..n-1 such that q[i] is the
	 *	  i-th vertex in a postorder numbering of the tree.
	 *
	 *        ( 2/7/95 modified by X.Li:
	 *          q is a vector indexed by 0:n-1 such that vertex i is the
	 *          q[i]-th vertex in a postorder numbering of the tree.
	 *          That is, this is the inverse of the previous q. )
	 *
	 *	In the child structure, lower-numbered children are represented
	 *	first, so that a tree which is already numbered in postorder
	 *	will not have its order changed.
	 *
	 *  Written by John Gilbert, Xerox, 10 Dec 1990.
	 *  Based on code written by John Gilbert at CMI in 1987.
	 */

	/**
	 * Depth-first search from vertex v.
	 */
	public static void etdfs (
		    int	  v,
		    int   first_kid[],
		    int   next_kid[],
		    int   post[],
		    int   postnum[]
		    )
	{
		int	w;

		for (w = first_kid[v]; w != -1; w = next_kid[w]) {
			etdfs (w, first_kid, next_kid, post, postnum);
		}
		/* post[postnum++] = v; in Matlab */
		post[v] = (postnum[0])++;    /* Modified by X. Li on 08/10/07 */
	}

	/**
	 * Depth-first search from vertex n.  No recursion.
	 * This routine was contributed by Cédric Doucet, CEDRAT Group, Meylan, France.
	 */
	public static void nr_etdfs (int n, int parent[],
		       int first_kid[], int next_kid[],
		       int post[], int postnum)
	{
	    int current = n, first, next;

	    while (postnum != n){

	        /* no kid for the current node */
	        first = first_kid[current];

	        /* no first kid for the current node */
	        if (first == -1){

	            /* numbering this node because it has no kid */
	            post[current] = postnum++;

	            /* looking for the next kid */
	            next = next_kid[current];

	            while (next == -1){

	                /* no more kids : back to the parent node */
	                current = parent[current];

	                /* numbering the parent node */
	                post[current] = postnum++;

	                /* get the next kid */
	                next = next_kid[current];
		    }

	            /* stopping criterion */
	            if (postnum==n+1) return;

	            /* updating current node */
	            current = next;
	        }
	        /* updating current node */
	        else {
	            current = first;
		}
	    }
	}

	/**
	 * Post order a tree
	 */
	@SuppressWarnings("unused")
	public static int[] TreePostorder(int n, int parent[]) {
        int	first_kid[], next_kid[];	/* Linked list of children.	*/
        int	post[], postnum[] = new int[1];
		int	v, dad;

		/* Allocate storage for working arrays and results	*/
		first_kid = 	mxCallocInt (n+1);
		next_kid  = 	mxCallocInt (n+1);
		post	  = 	mxCallocInt (n+1);

		/* Set up structure describing children */
		for (v = 0; v <= n; first_kid[v++] = -1);
		for (v = n-1; v >= 0; v--) {
			dad = parent[v];
			next_kid[v] = first_kid[dad];
			first_kid[dad] = v;
		}

		/* Depth-first search from dummy root vertex #n */
		postnum[0] = 0;
		if (false) {
		/* recursion */
		etdfs (n, first_kid, next_kid, post, postnum);
		} else {
		/* no recursion */
		nr_etdfs(n, parent, first_kid, next_kid, post, postnum[0]);
		}

		first_kid = null;
		next_kid = null;

		return post;
	}

}

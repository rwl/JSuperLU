package lu.jsuper;

public class Dlu_xerbla {

	/**  -- LAPACK auxiliary routine (version 2.0) --
	    Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
	    Courant Institute, Argonne National Lab, and Rice University
	    September 30, 1994


	 Purpose
	 =======

	 XERBLA  is an error handler for the LAPACK routines.
	 It is called by an LAPACK routine if an input parameter has an
	 invalid value.  A message is printed and execution stops.

	 Installers may consider modifying the STOP statement in order to
	 call system-specific exception-handling facilities.

	 Arguments
	 =========

	 SRNAME  (input) CHARACTER*6
	         The name of the routine which called XERBLA.

	 INFO    (input) INT
	         The position of the invalid parameter in the parameter list

	         of the calling routine.

	=====================================================================
	*/
	public static int xerbla_(String srname, int info) {
	    System.out.printf("** On entry to %6s, parameter number %2d had an illegal value\n",
	    		srname, info);
		return 0;
	}

}

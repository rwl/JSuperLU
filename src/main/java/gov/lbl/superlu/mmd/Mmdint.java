/*
 *  Produced by f2java.  f2java is part of the Fortran-
 *  -to-Java project at the University of Tennessee Netlib
 *  numerical software repository.
 *
 *  Original authorship for the BLAS and LAPACK numerical
 *  routines may be found in the Fortran source, available at
 *  http://www.netlib.org.
 *
 *  Fortran input file: genmmd.f
 *  f2java version: 0.8.1
 *
 */
package gov.lbl.superlu.mmd;

import java.lang.*;
import org.netlib.util.*;



public class Mmdint {

// C***************************************************************
// C***************************************************************
// C***     MMDINT ..... MULT MINIMUM DEGREE INITIALIZATION     ***
// C***************************************************************
// C***************************************************************
// C
// C     AUTHOR - JOSEPH W.H. LIU
// C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
// C
// C     PURPOSE - THIS ROUTINE PERFORMS INITIALIZATION FOR THE
// C        MULTIPLE ELIMINATION VERSION OF THE MINIMUM DEGREE
// C        ALGORITHM.
// C
// C     INPUT PARAMETERS -
// C        NEQNS  - NUMBER OF EQUATIONS.
// C        (XADJ,ADJNCY) - ADJACENCY STRUCTURE.
// C
// C     OUTPUT PARAMETERS -
// C        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
// C        QSIZE  - SIZE OF SUPERNODE (INITIALIZED TO ONE).
// C        LLIST  - LINKED LIST.
// C        MARKER - MARKER VECTOR.
// C
// C***************************************************************
// C
// C
// C***************************************************************
// C
// C         INTEGER*2  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
// C     1              LLIST(1) , MARKER(1), QSIZE(1)
// C
// C***************************************************************
// C

public static void mmdint (int neqns,
int [] xadj, int _xadj_offset,
int [] adjncy, int _adjncy_offset,
int [] dhead, int _dhead_offset,
int [] dforw, int _dforw_offset,
int [] dbakw, int _dbakw_offset,
int [] qsize, int _qsize_offset,
int [] llist, int _llist_offset,
int [] marker, int _marker_offset)  {

int fnode= 0;
int ndeg= 0;
int node= 0;
{
for (node = 1; node <= neqns; node++) {
dhead[(node-(1))+ _dhead_offset] = 0;
qsize[(node-(1))+ _qsize_offset] = 1;
marker[(node-(1))+ _marker_offset] = 0;
llist[(node-(1))+ _llist_offset] = 0;
Dummy.label("Mmdint",100);
}              //  Close for() loop.
}
// C        ------------------------------------------
// C        INITIALIZE THE DEGREE DOUBLY LINKED LISTS.
// C        ------------------------------------------
{
for (node = 1; node <= neqns; node++) {
ndeg = ((xadj[((node+1)-(1))+ _xadj_offset]-xadj[(node-(1))+ _xadj_offset])+1);
fnode = dhead[(ndeg-(1))+ _dhead_offset];
dforw[(node-(1))+ _dforw_offset] = fnode;
dhead[(ndeg-(1))+ _dhead_offset] = node;
if ((fnode > 0)) {
    dbakw[(fnode-(1))+ _dbakw_offset] = node;
}
    dbakw[(node-(1))+ _dbakw_offset] = (-(ndeg));
Dummy.label("Mmdint",200);
}              //  Close for() loop.
}
Dummy.go_to("Mmdint",999999);
// C
Dummy.label("Mmdint",999999);
return;
   }
} // End class.

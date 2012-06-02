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



public class Genmmd {

// C Sivan: I modified INTEGER*2 -> INTEGER*4
//
// C***************************************************************
// C***************************************************************
// C****     GENMMD ..... MULTIPLE MINIMUM EXTERNAL DEGREE     ****
// C***************************************************************
// C***************************************************************
// C
// C     AUTHOR - JOSEPH W.H. LIU
// C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
// C
// C     PURPOSE - THIS ROUTINE IMPLEMENTS THE MINIMUM DEGREE
// C        ALGORITHM.  IT MAKES USE OF THE IMPLICIT REPRESENTATION
// C        OF ELIMINATION GRAPHS BY QUOTIENT GRAPHS, AND THE
// C        NOTION OF INDISTINGUISHABLE NODES.  IT ALSO IMPLEMENTS
// C        THE MODIFICATIONS BY MULTIPLE ELIMINATION AND MINIMUM
// C        EXTERNAL DEGREE.
// C        ---------------------------------------------
// C        CAUTION - THE ADJACENCY VECTOR ADJNCY WILL BE
// C        DESTROYED.
// C        ---------------------------------------------
// C
// C     INPUT PARAMETERS -
// C        NEQNS  - NUMBER OF EQUATIONS.
// C        (XADJ,ADJNCY) - THE ADJACENCY STRUCTURE.
// C        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
// C        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) INTEGER
// C                 (ANY SMALLER ESTIMATE WILL DO) FOR MARKING
// C                 NODES.
// C
// C     OUTPUT PARAMETERS -
// C        PERM   - THE MINIMUM DEGREE ORDERING.
// C        INVP   - THE INVERSE OF PERM.
// C        NOFSUB - AN UPPER BOUND ON THE NUMBER OF NONZERO
// C                 SUBSCRIPTS FOR THE COMPRESSED STORAGE SCHEME.
// C
// C     WORKING PARAMETERS -
// C        DHEAD  - VECTOR FOR HEAD OF DEGREE LISTS.
// C        INVP   - USED TEMPORARILY FOR DEGREE FORWARD LINK.
// C        PERM   - USED TEMPORARILY FOR DEGREE BACKWARD LINK.
// C        QSIZE  - VECTOR FOR SIZE OF SUPERNODES.
// C        LLIST  - VECTOR FOR TEMPORARY LINKED LISTS.
// C        MARKER - A TEMPORARY MARKER VECTOR.
// C
// C     PROGRAM SUBROUTINES -
// C        MMDELM, MMDINT, MMDNUM, MMDUPD.
// C
// C***************************************************************
// C
// C
// C***************************************************************
// C
// C         INTEGER*2  ADJNCY(1), DHEAD(1) , INVP(1)  , LLIST(1) ,
// C     1              MARKER(1), PERM(1)  , QSIZE(1)
// C
// C***************************************************************
// C

public static void genmmd (int neqns,
int [] xadj, int _xadj_offset,
int [] adjncy, int _adjncy_offset,
int [] invp, int _invp_offset,
int [] perm, int _perm_offset,
int delta,
int [] dhead, int _dhead_offset,
int [] qsize, int _qsize_offset,
int [] llist, int _llist_offset,
int [] marker, int _marker_offset,
int maxint,
intW nofsub)  {

int ehead= 0;
int i= 0;
intW mdeg= new intW(0);
int mdlmt= 0;
int mdnode= 0;
int nextmd= 0;
int num= 0;
intW tag= new intW(0);
if ((neqns <= 0)) {
    Dummy.go_to("Genmmd",999999);
}
    // C
// C        ------------------------------------------------
// C        INITIALIZATION FOR THE MINIMUM DEGREE ALGORITHM.
// C        ------------------------------------------------
nofsub.val = 0;
Mmdint.mmdint(neqns,xadj,_xadj_offset,adjncy,_adjncy_offset,dhead,_dhead_offset,invp,_invp_offset,perm,_perm_offset,qsize,_qsize_offset,llist,_llist_offset,marker,_marker_offset);
// C
// C        ----------------------------------------------
// C        NUM COUNTS THE NUMBER OF ORDERED NODES PLUS 1.
// C        ----------------------------------------------
num = 1;
// C
// C        -----------------------------
// C        ELIMINATE ALL ISOLATED NODES.
// C        -----------------------------
nextmd = dhead[(1-(1))+ _dhead_offset];
label100:
   Dummy.label("Genmmd",100);
if ((nextmd <= 0)) {
    Dummy.go_to("Genmmd",200);
}
    mdnode = nextmd;
nextmd = invp[(mdnode-(1))+ _invp_offset];
marker[(mdnode-(1))+ _marker_offset] = maxint;
invp[(mdnode-(1))+ _invp_offset] = (-(num));
num = (num+1);
Dummy.go_to("Genmmd",100);
// C
label200:
   Dummy.label("Genmmd",200);
// C        ----------------------------------------
// C        SEARCH FOR NODE OF THE MINIMUM DEGREE.
// C        MDEG IS THE CURRENT MINIMUM DEGREE;
// C        TAG IS USED TO FACILITATE MARKING NODES.
// C        ----------------------------------------
if ((num > neqns)) {
    Dummy.go_to("Genmmd",1000);
}
    tag.val = 1;
dhead[(1-(1))+ _dhead_offset] = 0;
mdeg.val = 2;
label300:
   Dummy.label("Genmmd",300);
if ((dhead[(mdeg.val-(1))+ _dhead_offset] > 0)) {
    Dummy.go_to("Genmmd",400);
}
    mdeg.val = (mdeg.val+1);
Dummy.go_to("Genmmd",300);
label400:
   Dummy.label("Genmmd",400);
// C            -------------------------------------------------
// C            USE VALUE OF DELTA TO SET UP MDLMT, WHICH GOVERNS
// C            WHEN A DEGREE UPDATE IS TO BE PERFORMED.
// C            -------------------------------------------------
mdlmt = (mdeg.val+delta);
ehead = 0;
// C
label500:
   Dummy.label("Genmmd",500);
mdnode = dhead[(mdeg.val-(1))+ _dhead_offset];
if ((mdnode > 0)) {
    Dummy.go_to("Genmmd",600);
}
    mdeg.val = (mdeg.val+1);
if ((mdeg.val > mdlmt)) {
    Dummy.go_to("Genmmd",900);
}
    Dummy.go_to("Genmmd",500);
label600:
   Dummy.label("Genmmd",600);
// C                ----------------------------------------
// C                REMOVE MDNODE FROM THE DEGREE STRUCTURE.
// C                ----------------------------------------
nextmd = invp[(mdnode-(1))+ _invp_offset];
dhead[(mdeg.val-(1))+ _dhead_offset] = nextmd;
if ((nextmd > 0)) {
    perm[(nextmd-(1))+ _perm_offset] = (-(mdeg.val));
}
    invp[(mdnode-(1))+ _invp_offset] = (-(num));
nofsub.val = (((nofsub.val+mdeg.val)+qsize[(mdnode-(1))+ _qsize_offset])-2);
if (((num+qsize[(mdnode-(1))+ _qsize_offset]) > neqns)) {
    Dummy.go_to("Genmmd",1000);
}
    // C                ----------------------------------------------
// C                ELIMINATE MDNODE AND PERFORM QUOTIENT GRAPH
// C                TRANSFORMATION.  RESET TAG VALUE IF NECESSARY.
// C                ----------------------------------------------
tag.val = (tag.val+1);
if ((tag.val < maxint)) {
    Dummy.go_to("Genmmd",800);
}
    tag.val = 1;
{
for (i = 1; i <= neqns; i++) {
if ((marker[(i-(1))+ _marker_offset] < maxint)) {
    marker[(i-(1))+ _marker_offset] = 0;
}
    Dummy.label("Genmmd",700);
}              //  Close for() loop.
}
label800:
   Dummy.label("Genmmd",800);
Mmdelm.mmdelm(mdnode,xadj,_xadj_offset,adjncy,_adjncy_offset,dhead,_dhead_offset,invp,_invp_offset,perm,_perm_offset,qsize,_qsize_offset,llist,_llist_offset,marker,_marker_offset,maxint,tag.val);
num = (num+qsize[(mdnode-(1))+ _qsize_offset]);
llist[(mdnode-(1))+ _llist_offset] = ehead;
ehead = mdnode;
if ((delta >= 0)) {
    Dummy.go_to("Genmmd",500);
}
    label900:
   Dummy.label("Genmmd",900);
// C            -------------------------------------------
// C            UPDATE DEGREES OF THE NODES INVOLVED IN THE
// C            MINIMUM DEGREE NODES ELIMINATION.
// C            -------------------------------------------
if ((num > neqns)) {
    Dummy.go_to("Genmmd",1000);
}
    Mmdupd.mmdupd(ehead,neqns,xadj,_xadj_offset,adjncy,_adjncy_offset,delta,mdeg,dhead,_dhead_offset,invp,_invp_offset,perm,_perm_offset,qsize,_qsize_offset,llist,_llist_offset,marker,_marker_offset,maxint,tag);
Dummy.go_to("Genmmd",300);
// C
label1000:
   Dummy.label("Genmmd",1000);
Mmdnum.mmdnum(neqns,perm,_perm_offset,invp,_invp_offset,qsize,_qsize_offset);
Dummy.go_to("Genmmd",999999);
// C
Dummy.label("Genmmd",999999);
return;
   }
} // End class.

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



public class Mmdnum {

// C***************************************************************
// C***************************************************************
// C*****     MMDNUM ..... MULTI MINIMUM DEGREE NUMBERING     *****
// C***************************************************************
// C***************************************************************
// C
// C     AUTHOR - JOSEPH W.H. LIU
// C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
// C
// C     PURPOSE - THIS ROUTINE PERFORMS THE FINAL STEP IN
// C        PRODUCING THE PERMUTATION AND INVERSE PERMUTATION
// C        VECTORS IN THE MULTIPLE ELIMINATION VERSION OF THE
// C        MINIMUM DEGREE ORDERING ALGORITHM.
// C
// C     INPUT PARAMETERS -
// C        NEQNS  - NUMBER OF EQUATIONS.
// C        QSIZE  - SIZE OF SUPERNODES AT ELIMINATION.
// C
// C     UPDATED PARAMETERS -
// C        INVP   - INVERSE PERMUTATION VECTOR.  ON INPUT,
// C                 IF QSIZE(NODE)=0, THEN NODE HAS BEEN MERGED
// C                 INTO THE NODE -INVP(NODE); OTHERWISE,
// C                 -INVP(NODE) IS ITS INVERSE LABELLING.
// C
// C     OUTPUT PARAMETERS -
// C        PERM   - THE PERMUTATION VECTOR.
// C
// C***************************************************************
// C
// C
// C***************************************************************
// C
// C         INTEGER*2  INVP(1)  , PERM(1)  , QSIZE(1)
// C
// C***************************************************************
// C

public static void mmdnum (int neqns,
int [] perm, int _perm_offset,
int [] invp, int _invp_offset,
int [] qsize, int _qsize_offset)  {

int father= 0;
int nextf= 0;
int node= 0;
int nqsize= 0;
int num= 0;
int root= 0;
{
for (node = 1; node <= neqns; node++) {
nqsize = qsize[(node-(1))+ _qsize_offset];
if ((nqsize <= 0)) {
    perm[(node-(1))+ _perm_offset] = invp[(node-(1))+ _invp_offset];
}
    if ((nqsize > 0)) {
    perm[(node-(1))+ _perm_offset] = (-(invp[(node-(1))+ _invp_offset]));
}
    Dummy.label("Mmdnum",100);
}              //  Close for() loop.
}
// C        ------------------------------------------------------
// C        FOR EACH NODE WHICH HAS BEEN MERGED, DO THE FOLLOWING.
// C        ------------------------------------------------------
{
for (node = 1; node <= neqns; node++) {
if ((perm[(node-(1))+ _perm_offset] > 0)) {
    Dummy.go_to("Mmdnum",500);
}
    // C                -----------------------------------------
// C                TRACE THE MERGED TREE UNTIL ONE WHICH HAS
// C                NOT BEEN MERGED, CALL IT ROOT.
// C                -----------------------------------------
father = node;
label200:
   Dummy.label("Mmdnum",200);
if ((perm[(father-(1))+ _perm_offset] > 0)) {
    Dummy.go_to("Mmdnum",300);
}
    father = (-(perm[(father-(1))+ _perm_offset]));
Dummy.go_to("Mmdnum",200);
label300:
   Dummy.label("Mmdnum",300);
// C                -----------------------
// C                NUMBER NODE AFTER ROOT.
// C                -----------------------
root = father;
num = (perm[(root-(1))+ _perm_offset]+1);
invp[(node-(1))+ _invp_offset] = (-(num));
perm[(root-(1))+ _perm_offset] = num;
// C                ------------------------
// C                SHORTEN THE MERGED TREE.
// C                ------------------------
father = node;
label400:
   Dummy.label("Mmdnum",400);
nextf = (-(perm[(father-(1))+ _perm_offset]));
if ((nextf <= 0)) {
    Dummy.go_to("Mmdnum",500);
}
    perm[(father-(1))+ _perm_offset] = (-(root));
father = nextf;
Dummy.go_to("Mmdnum",400);
Dummy.label("Mmdnum",500);
}              //  Close for() loop.
}
// C        ----------------------
// C        READY TO COMPUTE PERM.
// C        ----------------------
{
for (node = 1; node <= neqns; node++) {
num = (-(invp[(node-(1))+ _invp_offset]));
invp[(node-(1))+ _invp_offset] = num;
perm[(num-(1))+ _perm_offset] = node;
Dummy.label("Mmdnum",600);
}              //  Close for() loop.
}
Dummy.go_to("Mmdnum",999999);
// C
Dummy.label("Mmdnum",999999);
return;
   }
} // End class.

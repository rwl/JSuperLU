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



public class Mmdupd {

// C***************************************************************
// C***************************************************************
// C*****     MMDUPD ..... MULTIPLE MINIMUM DEGREE UPDATE     *****
// C***************************************************************
// C***************************************************************
// C
// C     AUTHOR - JOSEPH W.H. LIU
// C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
// C
// C     PURPOSE - THIS ROUTINE UPDATES THE DEGREES OF NODES
// C        AFTER A MULTIPLE ELIMINATION STEP.
// C
// C     INPUT PARAMETERS -
// C        EHEAD  - THE BEGINNING OF THE LIST OF ELIMINATED
// C                 NODES (I.E., NEWLY FORMED ELEMENTS).
// C        NEQNS  - NUMBER OF EQUATIONS.
// C        (XADJ,ADJNCY) - ADJACENCY STRUCTURE.
// C        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
// C        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT)
// C                 INTEGER.
// C
// C     UPDATED PARAMETERS -
// C        MDEG   - NEW MINIMUM DEGREE AFTER DEGREE UPDATE.
// C        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
// C        QSIZE  - SIZE OF SUPERNODE.
// C        LLIST  - WORKING LINKED LIST.
// C        MARKER - MARKER VECTOR FOR DEGREE UPDATE.
// C        TAG    - TAG VALUE.
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

public static void mmdupd (int ehead,
int neqns,
int [] xadj, int _xadj_offset,
int [] adjncy, int _adjncy_offset,
int delta,
intW mdeg,
int [] dhead, int _dhead_offset,
int [] dforw, int _dforw_offset,
int [] dbakw, int _dbakw_offset,
int [] qsize, int _qsize_offset,
int [] llist, int _llist_offset,
int [] marker, int _marker_offset,
int maxint,
intW tag)  {

int deg= 0;
int deg0= 0;
int elmnt= 0;
int enode= 0;
int fnode= 0;
int i= 0;
int iq2= 0;
int istop= 0;
int istrt= 0;
int j= 0;
int jstop= 0;
int jstrt= 0;
int link= 0;
int mdeg0= 0;
int mtag= 0;
int nabor= 0;
int node= 0;
int q2head= 0;
int qxhead= 0;
mdeg0 = (mdeg.val+delta);
elmnt = ehead;
label100:
   Dummy.label("Mmdupd",100);
// C            -------------------------------------------------------
// C            FOR EACH OF THE NEWLY FORMED ELEMENT, DO THE FOLLOWING.
// C            (RESET TAG VALUE IF NECESSARY.)
// C            -------------------------------------------------------
if ((elmnt <= 0)) {
    Dummy.go_to("Mmdupd",999999);
}
    mtag = (tag.val+mdeg0);
if ((mtag < maxint)) {
    Dummy.go_to("Mmdupd",300);
}
    tag.val = 1;
{
for (i = 1; i <= neqns; i++) {
if ((marker[(i-(1))+ _marker_offset] < maxint)) {
    marker[(i-(1))+ _marker_offset] = 0;
}
    Dummy.label("Mmdupd",200);
}              //  Close for() loop.
}
mtag = (tag.val+mdeg0);
label300:
   Dummy.label("Mmdupd",300);
// C            ---------------------------------------------
// C            CREATE TWO LINKED LISTS FROM NODES ASSOCIATED
// C            WITH ELMNT: ONE WITH TWO NABORS (Q2HEAD) IN
// C            ADJACENCY STRUCTURE, AND THE OTHER WITH MORE
// C            THAN TWO NABORS (QXHEAD).  ALSO COMPUTE DEG0,
// C            NUMBER OF NODES IN THIS ELEMENT.
// C            ---------------------------------------------
q2head = 0;
qxhead = 0;
deg0 = 0;
link = elmnt;
label400:
   Dummy.label("Mmdupd",400);
istrt = xadj[(link-(1))+ _xadj_offset];
istop = (xadj[((link+1)-(1))+ _xadj_offset]-1);
{
for (i = istrt; i <= istop; i++) {
enode = adjncy[(i-(1))+ _adjncy_offset];
link = (-(enode));
{
  int _arif_tmp = enode;
if (_arif_tmp < 0)
      Dummy.go_to("Mmdupd",400);
else if (_arif_tmp == 0)
      Dummy.go_to("Mmdupd",800);
else   Dummy.go_to("Mmdupd",500);
}
// C
label500:
   Dummy.label("Mmdupd",500);
if ((qsize[(enode-(1))+ _qsize_offset] == 0)) {
    Dummy.go_to("Mmdupd",700);
}
    deg0 = (deg0+qsize[(enode-(1))+ _qsize_offset]);
marker[(enode-(1))+ _marker_offset] = mtag;
// C                        ----------------------------------
// C                        IF ENODE REQUIRES A DEGREE UPDATE,
// C                        THEN DO THE FOLLOWING.
// C                        ----------------------------------
if ((dbakw[(enode-(1))+ _dbakw_offset] != 0)) {
    Dummy.go_to("Mmdupd",700);
}
    // C                            ---------------------------------------
// C                            PLACE EITHER IN QXHEAD OR Q2HEAD LISTS.
// C                            ---------------------------------------
if ((dforw[(enode-(1))+ _dforw_offset] == 2)) {
    Dummy.go_to("Mmdupd",600);
}
    llist[(enode-(1))+ _llist_offset] = qxhead;
qxhead = enode;
Dummy.go_to("Mmdupd",700);
label600:
   Dummy.label("Mmdupd",600);
llist[(enode-(1))+ _llist_offset] = q2head;
q2head = enode;
Dummy.label("Mmdupd",700);
}              //  Close for() loop.
}
label800:
   Dummy.label("Mmdupd",800);
// C            --------------------------------------------
// C            FOR EACH ENODE IN Q2 LIST, DO THE FOLLOWING.
// C            --------------------------------------------
enode = q2head;
iq2 = 1;
label900:
   Dummy.label("Mmdupd",900);
if ((enode <= 0)) {
    Dummy.go_to("Mmdupd",1500);
}
    if ((dbakw[(enode-(1))+ _dbakw_offset] != 0)) {
    Dummy.go_to("Mmdupd",2200);
}
    tag.val = (tag.val+1);
deg = deg0;
// C                    ------------------------------------------
// C                    IDENTIFY THE OTHER ADJACENT ELEMENT NABOR.
// C                    ------------------------------------------
istrt = xadj[(enode-(1))+ _xadj_offset];
nabor = adjncy[(istrt-(1))+ _adjncy_offset];
if ((nabor == elmnt)) {
    nabor = adjncy[((istrt+1)-(1))+ _adjncy_offset];
}
    // C                    ------------------------------------------------
// C                    IF NABOR IS UNELIMINATED, INCREASE DEGREE COUNT.
// C                    ------------------------------------------------
link = nabor;
if ((dforw[(nabor-(1))+ _dforw_offset] < 0)) {
    Dummy.go_to("Mmdupd",1000);
}
    deg = (deg+qsize[(nabor-(1))+ _qsize_offset]);
Dummy.go_to("Mmdupd",2100);
label1000:
   Dummy.label("Mmdupd",1000);
// C                        --------------------------------------------
// C                        OTHERWISE, FOR EACH NODE IN THE 2ND ELEMENT,
// C                        DO THE FOLLOWING.
// C                        --------------------------------------------
istrt = xadj[(link-(1))+ _xadj_offset];
istop = (xadj[((link+1)-(1))+ _xadj_offset]-1);
{
for (i = istrt; i <= istop; i++) {
node = adjncy[(i-(1))+ _adjncy_offset];
link = (-(node));
if ((node == enode)) {
    Dummy.go_to("Mmdupd",1400);
}
    {
  int _arif_tmp = node;
if (_arif_tmp < 0)
      Dummy.go_to("Mmdupd",1000);
else if (_arif_tmp == 0)
      Dummy.go_to("Mmdupd",2100);
else   Dummy.go_to("Mmdupd",1100);
}
// C
label1100:
   Dummy.label("Mmdupd",1100);
if ((qsize[(node-(1))+ _qsize_offset] == 0)) {
    Dummy.go_to("Mmdupd",1400);
}
    if ((marker[(node-(1))+ _marker_offset] >= tag.val)) {
    Dummy.go_to("Mmdupd",1200);
}
    // C                                -------------------------------------
// C                                CASE WHEN NODE IS NOT YET CONSIDERED.
// C                                -------------------------------------
marker[(node-(1))+ _marker_offset] = tag.val;
deg = (deg+qsize[(node-(1))+ _qsize_offset]);
Dummy.go_to("Mmdupd",1400);
label1200:
   Dummy.label("Mmdupd",1200);
// C                            ----------------------------------------
// C                            CASE WHEN NODE IS INDISTINGUISHABLE FROM
// C                            ENODE.  MERGE THEM INTO A NEW SUPERNODE.
// C                            ----------------------------------------
if ((dbakw[(node-(1))+ _dbakw_offset] != 0)) {
    Dummy.go_to("Mmdupd",1400);
}
    if ((dforw[(node-(1))+ _dforw_offset] != 2)) {
    Dummy.go_to("Mmdupd",1300);
}
    qsize[(enode-(1))+ _qsize_offset] = (qsize[(enode-(1))+ _qsize_offset]+qsize[(node-(1))+ _qsize_offset]);
qsize[(node-(1))+ _qsize_offset] = 0;
marker[(node-(1))+ _marker_offset] = maxint;
dforw[(node-(1))+ _dforw_offset] = (-(enode));
dbakw[(node-(1))+ _dbakw_offset] = (-(maxint));
Dummy.go_to("Mmdupd",1400);
label1300:
   Dummy.label("Mmdupd",1300);
// C                            --------------------------------------
// C                            CASE WHEN NODE IS OUTMATCHED BY ENODE.
// C                            --------------------------------------
if ((dbakw[(node-(1))+ _dbakw_offset] == 0)) {
    dbakw[(node-(1))+ _dbakw_offset] = (-(maxint));
}
    Dummy.label("Mmdupd",1400);
}              //  Close for() loop.
}
Dummy.go_to("Mmdupd",2100);
label1500:
   Dummy.label("Mmdupd",1500);
// C                ------------------------------------------------
// C                FOR EACH ENODE IN THE QX LIST, DO THE FOLLOWING.
// C                ------------------------------------------------
enode = qxhead;
iq2 = 0;
label1600:
   Dummy.label("Mmdupd",1600);
if ((enode <= 0)) {
    Dummy.go_to("Mmdupd",2300);
}
    if ((dbakw[(enode-(1))+ _dbakw_offset] != 0)) {
    Dummy.go_to("Mmdupd",2200);
}
    tag.val = (tag.val+1);
deg = deg0;
// C                        ---------------------------------
// C                        FOR EACH UNMARKED NABOR OF ENODE,
// C                        DO THE FOLLOWING.
// C                        ---------------------------------
istrt = xadj[(enode-(1))+ _xadj_offset];
istop = (xadj[((enode+1)-(1))+ _xadj_offset]-1);
{
for (i = istrt; i <= istop; i++) {
nabor = adjncy[(i-(1))+ _adjncy_offset];
if ((nabor == 0)) {
    Dummy.go_to("Mmdupd",2100);
}
    if ((marker[(nabor-(1))+ _marker_offset] >= tag.val)) {
    Dummy.go_to("Mmdupd",2000);
}
    marker[(nabor-(1))+ _marker_offset] = tag.val;
link = nabor;
// C                                ------------------------------
// C                                IF UNELIMINATED, INCLUDE IT IN
// C                                DEG COUNT.
// C                                ------------------------------
if ((dforw[(nabor-(1))+ _dforw_offset] < 0)) {
    Dummy.go_to("Mmdupd",1700);
}
    deg = (deg+qsize[(nabor-(1))+ _qsize_offset]);
Dummy.go_to("Mmdupd",2000);
label1700:
   Dummy.label("Mmdupd",1700);
// C                                    -------------------------------
// C                                    IF ELIMINATED, INCLUDE UNMARKED
// C                                    NODES IN THIS ELEMENT INTO THE
// C                                    DEGREE COUNT.
// C                                    -------------------------------
jstrt = xadj[(link-(1))+ _xadj_offset];
jstop = (xadj[((link+1)-(1))+ _xadj_offset]-1);
{
for (j = jstrt; j <= jstop; j++) {
node = adjncy[(j-(1))+ _adjncy_offset];
link = (-(node));
{
  int _arif_tmp = node;
if (_arif_tmp < 0)
      Dummy.go_to("Mmdupd",1700);
else if (_arif_tmp == 0)
      Dummy.go_to("Mmdupd",2000);
else   Dummy.go_to("Mmdupd",1800);
}
// C
label1800:
   Dummy.label("Mmdupd",1800);
if ((marker[(node-(1))+ _marker_offset] >= tag.val)) {
    Dummy.go_to("Mmdupd",1900);
}
    marker[(node-(1))+ _marker_offset] = tag.val;
deg = (deg+qsize[(node-(1))+ _qsize_offset]);
Dummy.label("Mmdupd",1900);
}              //  Close for() loop.
}
Dummy.label("Mmdupd",2000);
}              //  Close for() loop.
}
label2100:
   Dummy.label("Mmdupd",2100);
// C                    -------------------------------------------
// C                    UPDATE EXTERNAL DEGREE OF ENODE IN DEGREE
// C                    STRUCTURE, AND MDEG (MIN DEG) IF NECESSARY.
// C                    -------------------------------------------
deg = ((deg-qsize[(enode-(1))+ _qsize_offset])+1);
fnode = dhead[(deg-(1))+ _dhead_offset];
dforw[(enode-(1))+ _dforw_offset] = fnode;
dbakw[(enode-(1))+ _dbakw_offset] = (-(deg));
if ((fnode > 0)) {
    dbakw[(fnode-(1))+ _dbakw_offset] = enode;
}
    dhead[(deg-(1))+ _dhead_offset] = enode;
if ((deg < mdeg.val)) {
    mdeg.val = deg;
}
    label2200:
   Dummy.label("Mmdupd",2200);
// C                    ----------------------------------
// C                    GET NEXT ENODE IN CURRENT ELEMENT.
// C                    ----------------------------------
enode = llist[(enode-(1))+ _llist_offset];
if ((iq2 == 1)) {
    Dummy.go_to("Mmdupd",900);
}
    Dummy.go_to("Mmdupd",1600);
label2300:
   Dummy.label("Mmdupd",2300);
// C            -----------------------------
// C            GET NEXT ELEMENT IN THE LIST.
// C            -----------------------------
tag.val = mtag;
elmnt = llist[(elmnt-(1))+ _llist_offset];
Dummy.go_to("Mmdupd",100);
// C
Dummy.label("Mmdupd",999999);
return;
   }
} // End class.

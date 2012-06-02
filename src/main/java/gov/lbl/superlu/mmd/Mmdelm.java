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



public class Mmdelm {

// C***************************************************************
// C***************************************************************
// C**     MMDELM ..... MULTIPLE MINIMUM DEGREE ELIMINATION     ***
// C***************************************************************
// C***************************************************************
// C
// C     AUTHOR - JOSEPH W.H. LIU
// C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
// C
// C     PURPOSE - THIS ROUTINE ELIMINATES THE NODE MDNODE OF
// C        MINIMUM DEGREE FROM THE ADJACENCY STRUCTURE, WHICH
// C        IS STORED IN THE QUOTIENT GRAPH FORMAT.  IT ALSO
// C        TRANSFORMS THE QUOTIENT GRAPH REPRESENTATION OF THE
// C        ELIMINATION GRAPH.
// C
// C     INPUT PARAMETERS -
// C        MDNODE - NODE OF MINIMUM DEGREE.
// C        MAXINT - ESTIMATE OF MAXIMUM REPRESENTABLE (SHORT)
// C                 INTEGER.
// C        TAG    - TAG VALUE.
// C
// C     UPDATED PARAMETERS -
// C        (XADJ,ADJNCY) - UPDATED ADJACENCY STRUCTURE.
// C        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
// C        QSIZE  - SIZE OF SUPERNODE.
// C        MARKER - MARKER VECTOR.
// C        LLIST  - TEMPORARY LINKED LIST OF ELIMINATED NABORS.
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
// C        -----------------------------------------------
// C        FIND REACHABLE SET AND PLACE IN DATA STRUCTURE.
// C        -----------------------------------------------

public static void mmdelm (int mdnode,
int [] xadj, int _xadj_offset,
int [] adjncy, int _adjncy_offset,
int [] dhead, int _dhead_offset,
int [] dforw, int _dforw_offset,
int [] dbakw, int _dbakw_offset,
int [] qsize, int _qsize_offset,
int [] llist, int _llist_offset,
int [] marker, int _marker_offset,
int maxint,
int tag)  {

int elmnt= 0;
int i= 0;
int istop= 0;
int istrt= 0;
int j= 0;
int jstop= 0;
int jstrt= 0;
int link= 0;
int nabor= 0;
int node= 0;
int npv= 0;
int nqnbrs= 0;
int nxnode= 0;
int pvnode= 0;
int rlmt= 0;
int rloc= 0;
int rnode= 0;
int xqnbr= 0;
marker[(mdnode-(1))+ _marker_offset] = tag;
istrt = xadj[(mdnode-(1))+ _xadj_offset];
istop = (xadj[((mdnode+1)-(1))+ _xadj_offset]-1);
// C        -------------------------------------------------------
// C        ELMNT POINTS TO THE BEGINNING OF THE LIST OF ELIMINATED
// C        NABORS OF MDNODE, AND RLOC GIVES THE STORAGE LOCATION
// C        FOR THE NEXT REACHABLE NODE.
// C        -------------------------------------------------------
elmnt = 0;
rloc = istrt;
rlmt = istop;
{
for (i = istrt; i <= istop; i++) {
nabor = adjncy[(i-(1))+ _adjncy_offset];
if ((nabor == 0)) {
    Dummy.go_to("Mmdelm",300);
}
    if ((marker[(nabor-(1))+ _marker_offset] >= tag)) {
    Dummy.go_to("Mmdelm",200);
}
    marker[(nabor-(1))+ _marker_offset] = tag;
if ((dforw[(nabor-(1))+ _dforw_offset] < 0)) {
    Dummy.go_to("Mmdelm",100);
}
    adjncy[(rloc-(1))+ _adjncy_offset] = nabor;
rloc = (rloc+1);
Dummy.go_to("Mmdelm",200);
label100:
   Dummy.label("Mmdelm",100);
llist[(nabor-(1))+ _llist_offset] = elmnt;
elmnt = nabor;
Dummy.label("Mmdelm",200);
}              //  Close for() loop.
}
label300:
   Dummy.label("Mmdelm",300);
// C            -----------------------------------------------------
// C            MERGE WITH REACHABLE NODES FROM GENERALIZED ELEMENTS.
// C            -----------------------------------------------------
if ((elmnt <= 0)) {
    Dummy.go_to("Mmdelm",1000);
}
    adjncy[(rlmt-(1))+ _adjncy_offset] = (-(elmnt));
link = elmnt;
label400:
   Dummy.label("Mmdelm",400);
jstrt = xadj[(link-(1))+ _xadj_offset];
jstop = (xadj[((link+1)-(1))+ _xadj_offset]-1);
{
for (j = jstrt; j <= jstop; j++) {
node = adjncy[(j-(1))+ _adjncy_offset];
link = (-(node));
{
  int _arif_tmp = node;
if (_arif_tmp < 0)
      Dummy.go_to("Mmdelm",400);
else if (_arif_tmp == 0)
      Dummy.go_to("Mmdelm",900);
else   Dummy.go_to("Mmdelm",500);
}
label500:
   Dummy.label("Mmdelm",500);
if (((marker[(node-(1))+ _marker_offset] >= tag) || (dforw[(node-(1))+ _dforw_offset] < 0))) {
    Dummy.go_to("Mmdelm",800);
}
    marker[(node-(1))+ _marker_offset] = tag;
// C                            ---------------------------------
// C                            USE STORAGE FROM ELIMINATED NODES
// C                            IF NECESSARY.
// C                            ---------------------------------
label600:
   Dummy.label("Mmdelm",600);
if ((rloc < rlmt)) {
    Dummy.go_to("Mmdelm",700);
}
    link = (-(adjncy[(rlmt-(1))+ _adjncy_offset]));
rloc = xadj[(link-(1))+ _xadj_offset];
rlmt = (xadj[((link+1)-(1))+ _xadj_offset]-1);
Dummy.go_to("Mmdelm",600);
label700:
   Dummy.label("Mmdelm",700);
adjncy[(rloc-(1))+ _adjncy_offset] = node;
rloc = (rloc+1);
Dummy.label("Mmdelm",800);
}              //  Close for() loop.
}
label900:
   Dummy.label("Mmdelm",900);
elmnt = llist[(elmnt-(1))+ _llist_offset];
Dummy.go_to("Mmdelm",300);
label1000:
   Dummy.label("Mmdelm",1000);
if ((rloc <= rlmt)) {
    adjncy[(rloc-(1))+ _adjncy_offset] = 0;
}
    // C        --------------------------------------------------------
// C        FOR EACH NODE IN THE REACHABLE SET, DO THE FOLLOWING ...
// C        --------------------------------------------------------
link = mdnode;
label1100:
   Dummy.label("Mmdelm",1100);
istrt = xadj[(link-(1))+ _xadj_offset];
istop = (xadj[((link+1)-(1))+ _xadj_offset]-1);
{
for (i = istrt; i <= istop; i++) {
rnode = adjncy[(i-(1))+ _adjncy_offset];
link = (-(rnode));
{
  int _arif_tmp = rnode;
if (_arif_tmp < 0)
      Dummy.go_to("Mmdelm",1100);
else if (_arif_tmp == 0)
      Dummy.go_to("Mmdelm",1800);
else   Dummy.go_to("Mmdelm",1200);
}
label1200:
   Dummy.label("Mmdelm",1200);
// C                --------------------------------------------
// C                IF RNODE IS IN THE DEGREE LIST STRUCTURE ...
// C                --------------------------------------------
pvnode = dbakw[(rnode-(1))+ _dbakw_offset];
if (((pvnode == 0) || (pvnode == ((-(maxint)))))) {
    Dummy.go_to("Mmdelm",1300);
}
    // C                    -------------------------------------
// C                    THEN REMOVE RNODE FROM THE STRUCTURE.
// C                    -------------------------------------
nxnode = dforw[(rnode-(1))+ _dforw_offset];
if ((nxnode > 0)) {
    dbakw[(nxnode-(1))+ _dbakw_offset] = pvnode;
}
    if ((pvnode > 0)) {
    dforw[(pvnode-(1))+ _dforw_offset] = nxnode;
}
    npv = (-(pvnode));
if ((pvnode < 0)) {
    dhead[(npv-(1))+ _dhead_offset] = nxnode;
}
    label1300:
   Dummy.label("Mmdelm",1300);
// C                ----------------------------------------
// C                PURGE INACTIVE QUOTIENT NABORS OF RNODE.
// C                ----------------------------------------
jstrt = xadj[(rnode-(1))+ _xadj_offset];
jstop = (xadj[((rnode+1)-(1))+ _xadj_offset]-1);
xqnbr = jstrt;
{
for (j = jstrt; j <= jstop; j++) {
nabor = adjncy[(j-(1))+ _adjncy_offset];
if ((nabor == 0)) {
    Dummy.go_to("Mmdelm",1500);
}
    if ((marker[(nabor-(1))+ _marker_offset] >= tag)) {
    Dummy.go_to("Mmdelm",1400);
}
    adjncy[(xqnbr-(1))+ _adjncy_offset] = nabor;
xqnbr = (xqnbr+1);
Dummy.label("Mmdelm",1400);
}              //  Close for() loop.
}
label1500:
   Dummy.label("Mmdelm",1500);
// C                ----------------------------------------
// C                IF NO ACTIVE NABOR AFTER THE PURGING ...
// C                ----------------------------------------
nqnbrs = (xqnbr-jstrt);
if ((nqnbrs > 0)) {
    Dummy.go_to("Mmdelm",1600);
}
    // C                    -----------------------------
// C                    THEN MERGE RNODE WITH MDNODE.
// C                    -----------------------------
qsize[(mdnode-(1))+ _qsize_offset] = (qsize[(mdnode-(1))+ _qsize_offset]+qsize[(rnode-(1))+ _qsize_offset]);
qsize[(rnode-(1))+ _qsize_offset] = 0;
marker[(rnode-(1))+ _marker_offset] = maxint;
dforw[(rnode-(1))+ _dforw_offset] = (-(mdnode));
dbakw[(rnode-(1))+ _dbakw_offset] = (-(maxint));
Dummy.go_to("Mmdelm",1700);
label1600:
   Dummy.label("Mmdelm",1600);
// C                --------------------------------------
// C                ELSE FLAG RNODE FOR DEGREE UPDATE, AND
// C                ADD MDNODE AS A NABOR OF RNODE.
// C                --------------------------------------
dforw[(rnode-(1))+ _dforw_offset] = (nqnbrs+1);
dbakw[(rnode-(1))+ _dbakw_offset] = 0;
adjncy[(xqnbr-(1))+ _adjncy_offset] = mdnode;
xqnbr = (xqnbr+1);
if ((xqnbr <= jstop)) {
    adjncy[(xqnbr-(1))+ _adjncy_offset] = 0;
}
    // C
Dummy.label("Mmdelm",1700);
}              //  Close for() loop.
}
label1800:
   Dummy.label("Mmdelm",1800);
Dummy.go_to("Mmdelm",999999);
// C
Dummy.label("Mmdelm",999999);
return;
   }
} // End class.

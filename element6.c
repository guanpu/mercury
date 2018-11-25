/* element6.f -- translated by f2c (version 20181026).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__6 = 6;
static integer c__9 = 9;
static integer c__250 = 250;
static integer c__3 = 3;
static integer c__80 = 80;
static integer c__23 = 23;
static integer c__5 = 5;
static doublereal c_b104 = 0.;
static doublereal c_b105 = 11239424.;
static doublereal c_b118 = 1.;
static doublereal c_b139 = 3.141592653589793;
static doublereal c_b141 = 6.2831853071795862;
static doublereal c_b149 = 360.;
static doublereal c_b198 = 10.;
static integer c__8 = 8;
static integer c__7 = 7;
static doublereal c_b245 = 1461.;
static doublereal c_b246 = 365.25;
static doublereal c_b299 = .3333333333333333;
static doublereal c_b317 = .33333333333333331;

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      ELEMENT6.FOR    (ErikSoft   5 June 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Makes output files containing Keplerian orbital elements from data created */
/* by Mercury6 and higher. */

/* The user specifies the names of the required objects in the file elements.in */
/* See subroutine M_FORMAT for the identities of each element in the EL array */
/* e.g. el(1)=a, el(2)=e etc. */

/* ------------------------------------------------------------------------------ */

/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_213[] = "(1x,a8,1x,f8.4,1x,f7.5,1x,f7.3,1p,e11.4,0p,1x,f"
	    "6.3,1x,f6.2)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;
    icilist ici__1;
    olist o__1;
    cllist cl__1;
    alist al__1, al__2;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer f_open(olist *), s_rsfe(cilist *), e_rsfe(void), f_clos(cllist *),
	     s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), f_back(alist *);
    double d_lg10(doublereal *);
    integer s_wsfi(icilist *), e_wsfi(void);
    double sqrt(doublereal), acos(doublereal), atan2(doublereal, doublereal), 
	    d_mod(doublereal *, doublereal *), d_sign(doublereal *, 
	    doublereal *), sin(doublereal), cos(doublereal);
    integer f_rew(alist *);

    /* Local variables */
    extern /* Subroutine */ int mco_ov2x__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), mco_iden__(doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *), mio_jd_y__(doublereal *, integer *, 
	    integer *, doublereal *), mce_spin__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), m_format__(char *, 
	    integer *, integer *, integer *, char *, char *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer line_num__;
    static doublereal a[2000];
    extern /* Subroutine */ int mxx_sort__(integer *, doublereal *, integer *)
	    ;
    static char c__[80*2000];
    static integer i__, j, k, l;
    static doublereal m[2000], s[3], v[6000]	/* was [3][2000] */, x[6000]	
	    /* was [3][2000] */;
    static char c1[1], c2[2], master_id__[8*2000];
    static integer firstflag;
    static doublereal t0, t1;
    static integer precision, timestyle;
    static char cc[80], id[8*2000];
    static doublereal tprevious, el[44000]	/* was [22][2000] */, gm, fr, 
	    is[2000], fv, vh[6000]	/* was [3][2000] */, xh[6000]	/* 
	    was [3][2000] */, ns[2000];
    static integer iel[22];
    static char fin[5];
    static integer nel;
    static char mem[80*200];
    static doublereal phi;
    static integer lim[200]	/* was [2][100] */, master_unit__[2000], code[
	    2000];
    static doublereal rfac;
    static integer nbig;
    static doublereal jcen[3];
    static integer nbod;
    static doublereal mcen, rcen;
    static integer lmem[200];
    static doublereal time;
    static integer year;
    static doublereal temp, vphi;
    static integer nsub;
    static doublereal rmax;
    static integer itmp, nsml;
    static char fout[250];
    static integer unit[2000];
    static logical test;
    static char type__[1];
    static integer nbig1, nbod1, iback[2000];
    static char check[1];
    static integer nchar, algor, lenin;
    static doublereal theta, teval;
    static integer nopen, nwait, month;
    static char style[1], header[250], infile__[250*50];
    static integer centre;
    static doublereal rhocgs, vtheta;
    static char string[250];
    extern /* Subroutine */ int mco_h2b__(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), mco_h2j__(doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static integer lenhead, allflag;
    extern /* Subroutine */ int mio_aei__(char *, char *, integer *, char *, 
	    integer *, char *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer ninfile;
    extern /* Subroutine */ int mio_err__(integer *, char *, integer *, char *
	    , integer *, char *, integer *, char *, integer *, ftnlen, ftnlen,
	     ftnlen, ftnlen), mio_spl__(integer *, char *, integer *, integer 
	    *, ftnlen);
    static integer nmaster;
    extern /* Subroutine */ int mco_h2cb__(doublereal *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern doublereal mio_c2fl__(char *, ftnlen), mio_c2re__(char *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    extern /* Subroutine */ int mco_x2el__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___6 = { 0, 14, 1, "(i3,1x,i2,1x,a80)", 0 };
    static cilist io___10 = { 0, 10, 0, "(a250)", 0 };
    static cilist io___15 = { 0, 10, 0, "(a250)", 0 };
    static cilist io___18 = { 0, 10, 0, "(a250)", 0 };
    static cilist io___21 = { 0, 10, 0, "(a250)", 0 };
    static cilist io___32 = { 0, 10, 1, "(a250)", 0 };
    static cilist io___39 = { 1, 10, 1, "(3a1)", 0 };
    static cilist io___43 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___44 = { 0, 10, 0, "(3x,i2,a62,i1)", 0 };
    static cilist io___56 = { 1, 10, 0, "(a)", 0 };
    static icilist io___61 = { 0, fin+2, 0, "(i2)", 2, 1 };
    static cilist io___71 = { 1, 10, 0, "(3x,a14)", 0 };
    static cilist io___74 = { 1, 10, 0, fin, 0 };
    static cilist io___77 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___94 = { 0, 0, 0, fout, 0 };
    static cilist io___95 = { 0, 0, 0, fout, 0 };
    static cilist io___96 = { 0, 6, 0, "(2a,/,a,i10)", 0 };
    static cilist io___97 = { 0, 10, 1, "(a1)", 0 };
    static cilist io___98 = { 0, 10, 0, "(/,a,f18.5,/)", 0 };
    static cilist io___99 = { 0, 10, 0, "(/,a,i10,1x,i2,1x,f8.5,/)", 0 };
    static cilist io___100 = { 0, 10, 0, "(/,a,f18.7,/)", 0 };
    static cilist io___101 = { 0, 10, 0, "(2a,/)", 0 };
    static cilist io___104 = { 0, 10, 0, fmt_213, 0 };



/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MERCURY.INC    (ErikSoft   4 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Parameters that you may want to alter at some point: */

/* NMAX  = maximum number of bodies */
/* CMAX  = maximum number of close-encounter minima monitored simultaneously */
/* NMESS = maximum number of messages in message.in */
/* HUGE  = an implausibly large number */
/* NFILES = maximum number of files that can be open at the same time */



/* ------------------------------------------------------------------------------ */

/* Constants: */

/* DR = conversion factor from degrees to radians */
/* K2 = Gaussian gravitational constant squared */
/* AU = astronomical unit in cm */
/* MSUN = mass of the Sun in g */



/* ------------------------------------------------------------------------------ */

    allflag = 0;
    tprevious = 0.;
    rhocgs = 498.06095345055081;

/* Read in output messages */
    ioin__1.inerr = 0;
    ioin__1.infilen = 10;
    ioin__1.infile = "message.in";
    ioin__1.inex = &test;
    ioin__1.inopen = 0;
    ioin__1.innum = 0;
    ioin__1.innamed = 0;
    ioin__1.inname = 0;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    if (! test) {
	s_wsfe(&io___5);
	do_fio(&c__1, " ERROR: This file is needed to continue: ", (ftnlen)41)
		;
	do_fio(&c__1, " message.in", (ftnlen)11);
	e_wsfe();
	s_stop("", (ftnlen)0);
    }
    o__1.oerr = 0;
    o__1.ounit = 14;
    o__1.ofnmlen = 10;
    o__1.ofnm = "message.in";
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
L10:
    i__1 = s_rsfe(&io___6);
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = do_fio(&c__1, (char *)&lmem[j - 1], (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = do_fio(&c__1, mem + (j - 1) * 80, (ftnlen)80);
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L20;
    }
    goto L10;
L20:
    cl__1.cerr = 0;
    cl__1.cunit = 14;
    cl__1.csta = 0;
    f_clos(&cl__1);

/* Open file containing parameters for this programme */
    ioin__1.inerr = 0;
    ioin__1.infilen = 10;
    ioin__1.infile = "element.in";
    ioin__1.inex = &test;
    ioin__1.inopen = 0;
    ioin__1.innum = 0;
    ioin__1.innamed = 0;
    ioin__1.inname = 0;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    if (test) {
	o__1.oerr = 0;
	o__1.ounit = 10;
	o__1.ofnmlen = 10;
	o__1.ofnm = "element.in";
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    } else {
	mio_err__(&c__6, mem + 6400, &lmem[80], mem + 6960, &lmem[87], " ", &
		c__1, "element.in", &c__9, (ftnlen)80, (ftnlen)80, (ftnlen)1, 
		(ftnlen)10);
    }

/* Read number of input files */
L30:
    s_rsfe(&io___10);
    do_fio(&c__1, string, (ftnlen)250);
    e_rsfe();
    if (*(unsigned char *)string == ')') {
	goto L30;
    }
    mio_spl__(&c__250, string, &nsub, lim, (ftnlen)250);
    i__1 = lim[(nsub << 1) - 2] - 1;
    ici__1.icierr = 0;
    ici__1.iciend = 0;
    ici__1.icirnum = 1;
    ici__1.icirlen = lim[(nsub << 1) - 1] - i__1;
    ici__1.iciunit = string + i__1;
    ici__1.icifmt = 0;
    s_rsli(&ici__1);
    do_lio(&c__3, &c__1, (char *)&ninfile, (ftnlen)sizeof(integer));
    e_rsli();

/* Make sure all the input files exist */
    i__1 = ninfile;
    for (j = 1; j <= i__1; ++j) {
L40:
	s_rsfe(&io___15);
	do_fio(&c__1, string, (ftnlen)250);
	e_rsfe();
	if (*(unsigned char *)string == ')') {
	    goto L40;
	}
	mio_spl__(&c__250, string, &nsub, lim, (ftnlen)250);
	i__2 = lim[0] - 1;
	s_copy(infile__ + (j - 1) * 250, string + i__2, lim[1] - lim[0] + 1, 
		lim[1] - i__2);
	ioin__1.inerr = 0;
	ioin__1.infilen = 250;
	ioin__1.infile = infile__ + (j - 1) * 250;
	ioin__1.inex = &test;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! test) {
	    mio_err__(&c__6, mem + 6400, &lmem[80], mem + 6960, &lmem[87], 
		    " ", &c__1, infile__ + (j - 1) * 250, &c__80, (ftnlen)80, 
		    (ftnlen)80, (ftnlen)1, (ftnlen)250);
	}
    }

/* What type elements does the user want? */
    centre = 0;
L45:
    s_rsfe(&io___18);
    do_fio(&c__1, string, (ftnlen)250);
    e_rsfe();
    if (*(unsigned char *)string == ')') {
	goto L45;
    }
    mio_spl__(&c__250, string, &nsub, lim, (ftnlen)250);
    i__1 = lim[(nsub << 1) - 2] - 1;
    s_copy(c2, string + i__1, (ftnlen)2, lim[(nsub << 1) - 2] + 1 - i__1);
    if (s_cmp(c2, "ce", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2, "CE", (ftnlen)
	    2, (ftnlen)2) == 0 || s_cmp(c2, "Ce", (ftnlen)2, (ftnlen)2) == 0) 
	    {
	centre = 0;
    } else if (s_cmp(c2, "ba", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2, "BA", (
	    ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2, "Ba", (ftnlen)2, (ftnlen)2)
	     == 0) {
	centre = 1;
    } else if (s_cmp(c2, "ja", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2, "JA", (
	    ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2, "Ja", (ftnlen)2, (ftnlen)2)
	     == 0) {
	centre = 2;
    } else {
	mio_err__(&c__6, mem + 6400, &lmem[80], mem + 8480, &lmem[106], " ", &
		c__1, "       Check element.in", &c__23, (ftnlen)80, (ftnlen)
		80, (ftnlen)1, (ftnlen)23);
    }

/* Read parameters used by this programme */
    timestyle = 1;
    for (j = 1; j <= 4; ++j) {
L50:
	s_rsfe(&io___21);
	do_fio(&c__1, string, (ftnlen)250);
	e_rsfe();
	if (*(unsigned char *)string == ')') {
	    goto L50;
	}
	mio_spl__(&c__250, string, &nsub, lim, (ftnlen)250);
	i__1 = lim[(nsub << 1) - 2] - 1;
	s_copy(c1, string + i__1, (ftnlen)1, lim[(nsub << 1) - 1] - i__1);
	if (j == 1) {
	    i__1 = lim[(nsub << 1) - 2] - 1;
	    ici__1.icierr = 0;
	    ici__1.iciend = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = lim[(nsub << 1) - 1] - i__1;
	    ici__1.iciunit = string + i__1;
	    ici__1.icifmt = 0;
	    s_rsli(&ici__1);
	    do_lio(&c__5, &c__1, (char *)&teval, (ftnlen)sizeof(doublereal));
	    e_rsli();
	}
	teval = abs(teval) * .999;
	if (j == 2 && (*(unsigned char *)c1 == 'd' || *(unsigned char *)c1 == 
		'D')) {
	    timestyle = 0;
	}
	if (j == 3 && (*(unsigned char *)c1 == 'y' || *(unsigned char *)c1 == 
		'Y')) {
	    timestyle += 2;
	}
	if (j == 4) {
	    m_format__(string, &timestyle, &nel, iel, fout, header, &lenhead, 
		    (ftnlen)250, (ftnlen)250, (ftnlen)250);
	}
    }

/* Read in the names of the objects for which orbital elements are required */
    nopen = 0;
    nwait = 0;
    nmaster = 0;
L60:
    i__1 = s_rsfe(&io___32);
    if (i__1 != 0) {
	goto L70;
    }
    i__1 = do_fio(&c__1, string, (ftnlen)250);
    if (i__1 != 0) {
	goto L70;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L70;
    }
    mio_spl__(&c__250, string, &nsub, lim, (ftnlen)250);
    if (*(unsigned char *)string == ')' || lim[0] == -1) {
	goto L60;
    }

/* Either open an aei file for this object or put it on the waiting list */
    ++nmaster;
/* Computing MIN */
    i__1 = 7, i__2 = lim[1] - lim[0];
    itmp = min(i__1,i__2);
    s_copy(master_id__ + (nmaster - 1 << 3), "        ", (ftnlen)8, (ftnlen)8)
	    ;
    i__1 = lim[0] - 1;
    s_copy(master_id__ + (nmaster - 1 << 3), string + i__1, itmp + 1, lim[0] 
	    + itmp - i__1);
    if (nopen < 50) {
	++nopen;
	master_unit__[nmaster - 1] = nopen + 10;
	mio_aei__(master_id__ + (nmaster - 1 << 3), ".aei", &master_unit__[
		nmaster - 1], header, &lenhead, mem, lmem, (ftnlen)8, (ftnlen)
		4, (ftnlen)250, (ftnlen)80);
    } else {
	++nwait;
	master_unit__[nmaster - 1] = -2;
    }
    goto L60;

L70:
/* If no objects are listed in ELEMENT.IN assume that all objects are required */
    if (nopen == 0) {
	allflag = 1;
    }
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);

/* ------------------------------------------------------------------------------ */

/*  LOOP  OVER  EACH  INPUT  FILE  CONTAINING  INTEGRATION  DATA */

L90:
    firstflag = 0;
    i__1 = ninfile;
    for (i__ = 1; i__ <= i__1; ++i__) {
	line_num__ = 0;
	o__1.oerr = 0;
	o__1.ounit = 10;
	o__1.ofnmlen = 250;
	o__1.ofnm = infile__ + (i__ - 1) * 250;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);

/* Loop over each time slice */
L100:
	++line_num__;
	i__2 = s_rsfe(&io___39);
	if (i__2 != 0) {
	    goto L100001;
	}
	i__2 = do_fio(&c__1, check, (ftnlen)1);
	if (i__2 != 0) {
	    goto L100001;
	}
	i__2 = do_fio(&c__1, style, (ftnlen)1);
	if (i__2 != 0) {
	    goto L100001;
	}
	i__2 = do_fio(&c__1, type__, (ftnlen)1);
	if (i__2 != 0) {
	    goto L100001;
	}
	i__2 = e_rsfe();
L100001:
	if (i__2 < 0) {
	    goto L900;
	}
	if (i__2 > 0) {
	    goto L666;
	}
	--line_num__;
	al__1.aerr = 0;
	al__1.aunit = 10;
	f_back(&al__1);

/* Check if this is an old style input file */
	if (*(unsigned char *)check == 12 && (*(unsigned char *)style == '0' 
		|| *(unsigned char *)style == '1' || *(unsigned char *)style 
		== '2' || *(unsigned char *)style == '3' || *(unsigned char *)
		style == '4')) {
	    s_wsfe(&io___43);
	    do_fio(&c__1, " ERROR: This is an old style data file", (ftnlen)
		    38);
	    do_fio(&c__1, "        Try running m_elem5.for instead.", (ftnlen)
		    40);
	    e_wsfe();
	    s_stop("", (ftnlen)0);
	}
	if (*(unsigned char *)check != 12) {
	    goto L666;
	}

/* ------------------------------------------------------------------------------ */

/*  IF  SPECIAL  INPUT,  READ  TIME,  PARAMETERS,  NAMES,  MASSES  ETC. */

	if (*(unsigned char *)type__ == 'a') {
	    ++line_num__;
	    s_rsfe(&io___44);
	    do_fio(&c__1, (char *)&algor, (ftnlen)sizeof(integer));
	    do_fio(&c__1, cc, (ftnlen)62);
	    do_fio(&c__1, (char *)&precision, (ftnlen)sizeof(integer));
	    e_rsfe();

/* Decompress the time, number of objects, central mass and J components etc. */
	    time = mio_c2fl__(cc, (ftnlen)8);
	    nbig = (integer) (mio_c2re__(cc + 8, &c_b104, &c_b105, &c__3, (
		    ftnlen)8) + .5);
	    nsml = (integer) (mio_c2re__(cc + 11, &c_b104, &c_b105, &c__3, (
		    ftnlen)8) + .5);
	    mcen = mio_c2fl__(cc + 14, (ftnlen)8);
	    jcen[0] = mio_c2fl__(cc + 22, (ftnlen)8);
	    jcen[1] = mio_c2fl__(cc + 30, (ftnlen)8);
	    jcen[2] = mio_c2fl__(cc + 38, (ftnlen)8);
	    rcen = mio_c2fl__(cc + 46, (ftnlen)8);
	    rmax = mio_c2fl__(cc + 54, (ftnlen)8);
	    d__1 = rmax / rcen;
	    rfac = d_lg10(&d__1);

/* Read in strings containing compressed data for each object */
	    i__2 = nbig + nsml;
	    for (j = 1; j <= i__2; ++j) {
		++line_num__;
		i__3 = s_rsfe(&io___56);
		if (i__3 != 0) {
		    goto L666;
		}
		i__3 = do_fio(&c__1, c__ + (j - 1) * 80, (ftnlen)51);
		if (i__3 != 0) {
		    goto L666;
		}
		i__3 = e_rsfe();
		if (i__3 != 0) {
		    goto L666;
		}
	    }

/* Create input format list */
	    if (precision == 1) {
		nchar = 2;
	    }
	    if (precision == 2) {
		nchar = 4;
	    }
	    if (precision == 3) {
		nchar = 7;
	    }
	    lenin = nchar * 6 + 3;
	    s_copy(fin, "(a00)", (ftnlen)5, (ftnlen)5);
	    s_wsfi(&io___61);
	    do_fio(&c__1, (char *)&lenin, (ftnlen)sizeof(integer));
	    e_wsfi();

/* For each object decompress its name, code number, mass, spin and density */
	    i__2 = nbig + nsml;
	    for (j = 1; j <= i__2; ++j) {
		k = (integer) (mio_c2re__(c__ + (j - 1) * 80, &c_b104, &
			c_b105, &c__3, (ftnlen)8) + .5);
		s_copy(id + (k - 1 << 3), c__ + ((j - 1) * 80 + 3), (ftnlen)8,
			 (ftnlen)8);
		el[k * 22 - 5] = mio_c2fl__(c__ + ((j - 1) * 80 + 11), (
			ftnlen)8);
		s[0] = mio_c2fl__(c__ + ((j - 1) * 80 + 19), (ftnlen)8);
		s[1] = mio_c2fl__(c__ + ((j - 1) * 80 + 27), (ftnlen)8);
		s[2] = mio_c2fl__(c__ + ((j - 1) * 80 + 35), (ftnlen)8);
		el[k * 22 - 2] = mio_c2fl__(c__ + ((j - 1) * 80 + 43), (
			ftnlen)8);

/* Calculate spin rate and longitude & inclination of spin vector */
		temp = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
		if (temp > 0.) {
		    d__1 = el[k * 22 - 5] * 2.959122082855911e-4;
		    d__2 = temp * 2.959122082855911e-4;
		    d__3 = el[k * 22 - 2] * rhocgs;
		    mce_spin__(&c_b118, &d__1, &d__2, &d__3, &el[k * 22 - 3]);
		    temp = s[2] / temp;
		    if (abs(temp) < 1.) {
			is[k - 1] = acos(temp);
			ns[k - 1] = atan2(s[0], -s[1]);
		    } else {
			if (temp > 0.) {
			    is[k - 1] = 0.;
			}
			if (temp < 0.) {
			    is[k - 1] = 3.141592653589793;
			}
			ns[k - 1] = 0.;
		    }
		} else {
		    el[k * 22 - 3] = 0.;
		    is[k - 1] = 0.;
		    ns[k - 1] = 0.;
		}

/* Find the object on the master list */
		unit[k - 1] = 0;
		i__3 = nmaster;
		for (l = 1; l <= i__3; ++l) {
		    if (s_cmp(id + (k - 1 << 3), master_id__ + (l - 1 << 3), (
			    ftnlen)8, (ftnlen)8) == 0) {
			unit[k - 1] = master_unit__[l - 1];
		    }
		}

/* If object is not on the master list, add it to the list now */
		if (unit[k - 1] == 0) {
		    ++nmaster;
		    s_copy(master_id__ + (nmaster - 1 << 3), id + (k - 1 << 3)
			    , (ftnlen)8, (ftnlen)8);

/* Either open an aei file for this object or put it on the waiting list */
		    if (allflag == 1) {
			if (nopen < 50) {
			    ++nopen;
			    master_unit__[nmaster - 1] = nopen + 10;
			    mio_aei__(master_id__ + (nmaster - 1 << 3), ".aei"
				    , &master_unit__[nmaster - 1], header, &
				    lenhead, mem, lmem, (ftnlen)8, (ftnlen)4, 
				    (ftnlen)250, (ftnlen)80);
			} else {
			    ++nwait;
			    master_unit__[nmaster - 1] = -2;
			}
		    } else {
			master_unit__[nmaster - 1] = -1;
		    }
		    unit[k - 1] = master_unit__[nmaster - 1];
		}
	    }

/* ------------------------------------------------------------------------------ */

/*  IF  NORMAL  INPUT,  READ  COMPRESSED  ORBITAL  VARIABLES  FOR  ALL  OBJECTS */

	} else if (*(unsigned char *)type__ == 'b') {
	    ++line_num__;
	    i__2 = s_rsfe(&io___71);
	    if (i__2 != 0) {
		goto L666;
	    }
	    i__2 = do_fio(&c__1, cc, (ftnlen)14);
	    if (i__2 != 0) {
		goto L666;
	    }
	    i__2 = e_rsfe();
	    if (i__2 != 0) {
		goto L666;
	    }

/* Decompress the time and the number of objects */
	    time = mio_c2fl__(cc, (ftnlen)8);
	    nbig = (integer) (mio_c2re__(cc + 8, &c_b104, &c_b105, &c__3, (
		    ftnlen)8) + .5);
	    nsml = (integer) (mio_c2re__(cc + 11, &c_b104, &c_b105, &c__3, (
		    ftnlen)8) + .5);
	    nbod = nbig + nsml;
	    if (firstflag == 0) {
		t0 = time;
	    }

/* Read in strings containing compressed data for each object */
	    i__2 = nbod;
	    for (j = 1; j <= i__2; ++j) {
		++line_num__;
		i__3 = s_rsfe(&io___74);
		if (i__3 != 0) {
		    goto L666;
		}
		i__3 = do_fio(&c__1, c__ + (j - 1) * 80, lenin);
		if (i__3 != 0) {
		    goto L666;
		}
		i__3 = e_rsfe();
		if (i__3 != 0) {
		    goto L666;
		}
	    }

/* Look for objects for which orbital elements are required */
	    m[0] = mcen * 2.959122082855911e-4;
	    i__2 = nbod;
	    for (j = 1; j <= i__2; ++j) {
		code[j - 1] = (integer) (mio_c2re__(c__ + (j - 1) * 80, &
			c_b104, &c_b105, &c__3, (ftnlen)8) + .5);
		if (code[j - 1] > 2000) {
		    s_wsfe(&io___77);
		    do_fio(&c__1, mem + 6400, lmem[80]);
		    do_fio(&c__1, mem + 7120, lmem[89]);
		    e_wsfe();
		    s_stop("", (ftnlen)0);
		}

/* Decompress orbital variables for each object */
		l = j + 1;
		m[l - 1] = el[code[j - 1] * 22 - 5] * 2.959122082855911e-4;
		fr = mio_c2re__(c__ + ((j - 1) * 80 + 3), &c_b104, &rfac, &
			nchar, (ftnlen)8);
		i__3 = nchar + 3;
		theta = mio_c2re__(c__ + ((j - 1) * 80 + i__3), &c_b104, &
			c_b139, &nchar, nchar + 11 - i__3);
		i__3 = (nchar << 1) + 3;
		phi = mio_c2re__(c__ + ((j - 1) * 80 + i__3), &c_b104, &
			c_b141, &nchar, (nchar << 1) + 11 - i__3);
		i__3 = nchar * 3 + 3;
		fv = mio_c2re__(c__ + ((j - 1) * 80 + i__3), &c_b104, &c_b118,
			 &nchar, nchar * 3 + 11 - i__3);
		i__3 = (nchar << 2) + 3;
		vtheta = mio_c2re__(c__ + ((j - 1) * 80 + i__3), &c_b104, &
			c_b139, &nchar, (nchar << 2) + 11 - i__3);
		i__3 = nchar * 5 + 3;
		vphi = mio_c2re__(c__ + ((j - 1) * 80 + i__3), &c_b104, &
			c_b141, &nchar, nchar * 5 + 11 - i__3);
		mco_ov2x__(&rcen, &rmax, m, &m[l - 1], &fr, &theta, &phi, &fv,
			 &vtheta, &vphi, &x[l * 3 - 3], &x[l * 3 - 2], &x[l * 
			3 - 1], &v[l * 3 - 3], &v[l * 3 - 2], &v[l * 3 - 1]);
		el[code[j - 1] * 22 - 7] = sqrt(x[l * 3 - 3] * x[l * 3 - 3] + 
			x[l * 3 - 2] * x[l * 3 - 2] + x[l * 3 - 1] * x[l * 3 
			- 1]);
	    }

/* Convert to barycentric, Jacobi or close-binary coordinates if desired */
	    nbod1 = nbod + 1;
	    nbig1 = nbig + 1;
	    mco_iden__(jcen, &nbod1, &nbig1, &temp, m, x, v, xh, vh);
	    if (centre == 1) {
		mco_h2b__(jcen, &nbod1, &nbig1, &temp, m, xh, vh, x, v);
	    }
	    if (centre == 2) {
		mco_h2j__(jcen, &nbod1, &nbig1, &temp, m, xh, vh, x, v);
	    }
	    if (centre == 0 && algor == 11) {
		mco_h2cb__(jcen, &nbod1, &nbig1, &temp, m, xh, vh, x, v);
	    }

/* Put Cartesian coordinates into element arrays */
	    i__2 = nbod;
	    for (j = 1; j <= i__2; ++j) {
		k = code[j - 1];
		l = j + 1;
		el[k * 22 - 13] = x[l * 3 - 3];
		el[k * 22 - 12] = x[l * 3 - 2];
		el[k * 22 - 11] = x[l * 3 - 1];
		el[k * 22 - 10] = v[l * 3 - 3];
		el[k * 22 - 9] = v[l * 3 - 2];
		el[k * 22 - 8] = v[l * 3 - 1];

/* Convert to Keplerian orbital elements */
		gm = (mcen + el[k * 22 - 5]) * 2.959122082855911e-4;
		mco_x2el__(&gm, &el[k * 22 - 13], &el[k * 22 - 12], &el[k * 
			22 - 11], &el[k * 22 - 10], &el[k * 22 - 9], &el[k * 
			22 - 8], &el[k * 22 - 15], &el[k * 22 - 21], &el[k * 
			22 - 20], &el[k * 22 - 16], &el[k * 22 - 18], &el[k * 
			22 - 17]);
		el[k * 22 - 22] = el[k * 22 - 15] / (1. - el[k * 22 - 21]);
		el[k * 22 - 14] = el[k * 22 - 22] * (el[k * 22 - 21] + 1.);
		d__1 = el[k * 22 - 16] - el[k * 22 - 18] + 6.2831853071795862;
		el[k * 22 - 19] = d_mod(&d__1, &c_b141);
/* Calculate true anomaly */
		if (el[k * 22 - 21] == 0.) {
		    el[k * 22 - 6] = el[k * 22 - 17];
		} else {
		    temp = (el[k * 22 - 15] * (el[k * 22 - 21] + 1.) / el[k * 
			    22 - 7] - 1.) / el[k * 22 - 21];
/* Computing MIN */
		    d__2 = abs(temp);
		    d__1 = min(d__2,1.);
		    temp = d_sign(&d__1, &temp);
		    el[k * 22 - 6] = acos(temp);
		    if (sin(el[k * 22 - 17]) < 0.) {
			el[k * 22 - 6] = 6.2831853071795862 - el[k * 22 - 6];
		    }
		}
/* Calculate obliquity */
		el[k * 22 - 4] = acos(cos(el[k * 22 - 20]) * cos(is[k - 1]) + 
			sin(el[k * 22 - 20]) * sin(is[k - 1]) * cos(ns[k - 1] 
			- el[k * 22 - 18]));

/* Convert angular elements from radians to degrees */
		for (l = 3; l <= 7; ++l) {
		    d__1 = el[l + k * 22 - 23] / .017453292519943295;
		    el[l + k * 22 - 23] = d_mod(&d__1, &c_b149);
		}
		el[k * 22 - 6] /= .017453292519943295;
		el[k * 22 - 4] /= .017453292519943295;
	    }

/* Convert time to desired format */
	    if (timestyle == 0) {
		t1 = time;
	    }
	    if (timestyle == 1) {
		mio_jd_y__(&time, &year, &month, &t1);
	    }
	    if (timestyle == 2) {
		t1 = time - t0;
	    }
	    if (timestyle == 3) {
		t1 = (time - t0) / 365.25;
	    }

/* If output is required at this epoch, write elements to appropriate files */
	    if (firstflag == 0 || (d__1 = time - tprevious, abs(d__1)) >= 
		    teval) {
		firstflag = 1;
		tprevious = time;

/* Write required elements to the appropriate aei file */
		i__2 = nbod;
		for (j = 1; j <= i__2; ++j) {
		    k = code[j - 1];
		    if (unit[k - 1] >= 10) {
			if (timestyle == 1) {
			    io___94.ciunit = unit[k - 1];
			    s_wsfe(&io___94);
			    do_fio(&c__1, (char *)&year, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&month, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(
				    doublereal));
			    i__3 = nel;
			    for (l = 1; l <= i__3; ++l) {
				do_fio(&c__1, (char *)&el[iel[l - 1] + k * 22 
					- 23], (ftnlen)sizeof(doublereal));
			    }
			    e_wsfe();
			} else {
			    io___95.ciunit = unit[k - 1];
			    s_wsfe(&io___95);
			    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(
				    doublereal));
			    i__3 = nel;
			    for (l = 1; l <= i__3; ++l) {
				do_fio(&c__1, (char *)&el[iel[l - 1] + k * 22 
					- 23], (ftnlen)sizeof(doublereal));
			    }
			    e_wsfe();
			}
		    }
		}
	    }

/* ------------------------------------------------------------------------------ */

/*  IF  TYPE  IS  NOT  'a'  OR  'b',  THE  INPUT  FILE  IS  CORRUPTED */

	} else {
	    goto L666;
	}

/* Move on to the next time slice */
	goto L100;

/* If input file is corrupted, try to continue from next uncorrupted time slice */
L666:
	s_wsfe(&io___96);
	do_fio(&c__1, mem + 9600, lmem[120]);
	do_fio(&c__1, infile__ + (i__ - 1) * 250, (ftnlen)60);
	do_fio(&c__1, mem + 8240, lmem[103]);
	do_fio(&c__1, (char *)&line_num__, (ftnlen)sizeof(integer));
	e_wsfe();
	*(unsigned char *)c1 = ' ';
	while(*(unsigned char *)c1 != 12) {
	    ++line_num__;
	    i__2 = s_rsfe(&io___97);
	    if (i__2 != 0) {
		goto L900;
	    }
	    i__2 = do_fio(&c__1, c1, (ftnlen)1);
	    if (i__2 != 0) {
		goto L900;
	    }
	    i__2 = e_rsfe();
	    if (i__2 != 0) {
		goto L900;
	    }
	}
	--line_num__;
	al__1.aerr = 0;
	al__1.aunit = 10;
	f_back(&al__1);

/* Move on to the next file containing integration data */
L900:
	cl__1.cerr = 0;
	cl__1.cunit = 10;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/* Close aei files */
    i__1 = nopen;
    for (j = 1; j <= i__1; ++j) {
	cl__1.cerr = 0;
	cl__1.cunit = j + 10;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    nopen = 0;

/* If some objects remain on waiting list, read through input files again */
    if (nwait > 0) {
	i__1 = nmaster;
	for (j = 1; j <= i__1; ++j) {
	    if (master_unit__[j - 1] >= 10) {
		master_unit__[j - 1] = -1;
	    }
	    if (master_unit__[j - 1] == -2 && nopen < 50) {
		++nopen;
		--nwait;
		master_unit__[j - 1] = nopen + 10;
		mio_aei__(master_id__ + (j - 1 << 3), ".aei", &master_unit__[
			j - 1], header, &lenhead, mem, lmem, (ftnlen)8, (
			ftnlen)4, (ftnlen)250, (ftnlen)80);
	    }
	}
	goto L90;
    }

/* ------------------------------------------------------------------------------ */

/*  CREATE  A  SUMMARY  OF  FINAL  MASSES  AND  ELEMENTS */

    o__1.oerr = 0;
    o__1.ounit = 10;
    o__1.ofnmlen = 11;
    o__1.ofnm = "element.out";
    o__1.orl = 0;
    o__1.osta = "unknown";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    al__2.aerr = 0;
    al__2.aunit = 10;
    f_rew(&al__2);

    if (timestyle == 0 || timestyle == 2) {
	s_wsfe(&io___98);
	do_fio(&c__1, " Time (days): ", (ftnlen)14);
	do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else if (timestyle == 1) {
	s_wsfe(&io___99);
	do_fio(&c__1, " Date: ", (ftnlen)7);
	do_fio(&c__1, (char *)&year, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&month, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else if (timestyle == 3) {
	s_wsfe(&io___100);
	do_fio(&c__1, " Time (years): ", (ftnlen)15);
	do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    s_wsfe(&io___101);
    do_fio(&c__1, "              a        e       i      mass", (ftnlen)42);
    do_fio(&c__1, "    Rot/day  Obl", (ftnlen)16);
    e_wsfe();

/* Sort surviving objects in order of increasing semi-major axis */
    i__1 = nbod;
    for (j = 1; j <= i__1; ++j) {
	k = code[j - 1];
	a[j - 1] = el[k * 22 - 22];
    }
    mxx_sort__(&nbod, a, iback);

/* Write values of a, e, i and m for surviving objects in an output file */
    i__1 = nbod;
    for (j = 1; j <= i__1; ++j) {
	k = code[iback[j - 1] - 1];
	s_wsfe(&io___104);
	do_fio(&c__1, id + (k - 1 << 3), (ftnlen)8);
	do_fio(&c__1, (char *)&el[k * 22 - 22], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&el[k * 22 - 21], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&el[k * 22 - 20], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&el[k * 22 - 5], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&el[k * 22 - 3], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&el[k * 22 - 4], (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/* ------------------------------------------------------------------------------ */

/* Format statements */

    return 0;
} /* MAIN__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_OV2X.FOR    (ErikSoft   28 February 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Converts output variables for an object to coordinates and velocities. */
/* The output variables are: */
/*  r = the radial distance */
/*  theta = polar angle */
/*  phi = azimuthal angle */
/*  fv = 1 / [1 + 2(ke/be)^2], where be and ke are the object's binding and */
/*                             kinetic energies. (Note that 0 < fv < 1). */
/*  vtheta = polar angle of velocity vector */
/*  vphi = azimuthal angle of the velocity vector */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_ov2x__(doublereal *rcen, doublereal *rmax, 
	doublereal *mcen, doublereal *m, doublereal *fr, doublereal *theta, 
	doublereal *phi, doublereal *fv, doublereal *vtheta, doublereal *vphi,
	 doublereal *x, doublereal *y, doublereal *z__, doublereal *u, 
	doublereal *v, doublereal *w)
{
    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), sin(
	    doublereal), cos(doublereal);

    /* Local variables */
    static doublereal r__, v1, temp;



/* Input/Output */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MERCURY.INC    (ErikSoft   4 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Parameters that you may want to alter at some point: */

/* NMAX  = maximum number of bodies */
/* CMAX  = maximum number of close-encounter minima monitored simultaneously */
/* NMESS = maximum number of messages in message.in */
/* HUGE  = an implausibly large number */
/* NFILES = maximum number of files that can be open at the same time */



/* ------------------------------------------------------------------------------ */

/* Constants: */

/* DR = conversion factor from degrees to radians */
/* K2 = Gaussian gravitational constant squared */
/* AU = astronomical unit in cm */
/* MSUN = mass of the Sun in g */



/* Local */

/* ------------------------------------------------------------------------------ */

    r__ = *rcen * pow_dd(&c_b198, fr);
    temp = sqrt((1. / *fv - 1.) * .5);
    v1 = sqrt(temp * 2. * (*mcen + *m) / r__);

    *x = r__ * sin(*theta) * cos(*phi);
    *y = r__ * sin(*theta) * sin(*phi);
    *z__ = r__ * cos(*theta);
    *u = v1 * sin(*vtheta) * cos(*vphi);
    *v = v1 * sin(*vtheta) * sin(*vphi);
    *w = v1 * cos(*vtheta);

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_ov2x__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCE_SPIN.FOR    (ErikSoft  2 December 1999) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates the spin rate (in rotations per day) for a fluid body given */
/* its mass, spin angular momentum and density. The routine assumes the */
/* body is a MacClaurin ellipsoid, whose axis ratio is defined by the */
/* quantity SS = SQRT(A^2/C^2 - 1), where A and C are the */
/* major and minor axes. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mce_spin__(doublereal *g, doublereal *mass, doublereal *
	spin, doublereal *rho, doublereal *rote)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal f;
    static integer k;
    static doublereal z__, s2, df, t23, dz, ss, tmp0, tmp1;
    extern /* Subroutine */ int m_sfunc__(doublereal *, doublereal *, 
	    doublereal *);



/* Input/Output */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MERCURY.INC    (ErikSoft   4 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Parameters that you may want to alter at some point: */

/* NMAX  = maximum number of bodies */
/* CMAX  = maximum number of close-encounter minima monitored simultaneously */
/* NMESS = maximum number of messages in message.in */
/* HUGE  = an implausibly large number */
/* NFILES = maximum number of files that can be open at the same time */



/* ------------------------------------------------------------------------------ */

/* Constants: */

/* DR = conversion factor from degrees to radians */
/* K2 = Gaussian gravitational constant squared */
/* AU = astronomical unit in cm */
/* MSUN = mass of the Sun in g */



/* Local */

/* ------------------------------------------------------------------------------ */

    t23 = .66666666666666663;
/* Computing 5th power */
    d__2 = *mass, d__3 = d__2, d__2 *= d__2;
    d__1 = *rho * 2467.4011002723391 * *rho / (d__3 * (d__2 * d__2) * 9.);
    tmp1 = *spin * *spin / (*rho * 6.2831853071795862 * *g) * pow_dd(&d__1, &
	    t23);

/* Calculate SS using Newton's method */
    ss = 1.;
    for (k = 1; k <= 20; ++k) {
	s2 = ss * ss;
	d__1 = s2 + 1.;
	tmp0 = pow_dd(&d__1, &t23);
	m_sfunc__(&ss, &z__, &dz);
	f = z__ * tmp0 - tmp1;
	df = tmp0 * (dz + ss * 4. * z__ / ((s2 + 1.) * 3.));
	ss -= f / df;
    }

    *rote = sqrt(*g * 6.2831853071795862 * *rho * z__) / 6.2831853071795862;

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mce_spin__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_EL2X.FOR    (ErikSoft  7 July 1999) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates Cartesian coordinates and velocities given Keplerian orbital */
/* elements (for elliptical, parabolic or hyperbolic orbits). */

/* Based on a routine from Levison and Duncan's SWIFT integrator. */

/*  mu = grav const * (central + secondary mass) */
/*  q = perihelion distance */
/*  e = eccentricity */
/*  i = inclination                 ) */
/*  p = longitude of perihelion !!! )   in */
/*  n = longitude of ascending node ) radians */
/*  l = mean anomaly                ) */

/*  x,y,z = Cartesian positions  ( units the same as a ) */
/*  u,v,w =     "     velocities ( units the same as sqrt(mu/a) ) */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_el2x__(doublereal *mu, doublereal *q, doublereal *e, 
	doublereal *i__, doublereal *p, doublereal *n, doublereal *l, 
	doublereal *x, doublereal *y, doublereal *z__, doublereal *u, 
	doublereal *v, doublereal *w)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int mco_sine__(doublereal *, doublereal *, 
	    doublereal *), mco_sinh__(doublereal *, doublereal *, doublereal *
	    );
    static doublereal a, g, z1, z2, z3, z4, d11, d12, ce, d13, cg, d21, ci, 
	    d22, d23, cn, se, sg, si, sn;
    extern doublereal orbel_zget__(doublereal *);
    static doublereal temp, romes;
    extern doublereal orbel_fhybrid__(doublereal *, doublereal *), mco_kep__(
	    doublereal *, doublereal *);



/* Input/Output */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MERCURY.INC    (ErikSoft   4 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Parameters that you may want to alter at some point: */

/* NMAX  = maximum number of bodies */
/* CMAX  = maximum number of close-encounter minima monitored simultaneously */
/* NMESS = maximum number of messages in message.in */
/* HUGE  = an implausibly large number */
/* NFILES = maximum number of files that can be open at the same time */



/* ------------------------------------------------------------------------------ */

/* Constants: */

/* DR = conversion factor from degrees to radians */
/* K2 = Gaussian gravitational constant squared */
/* AU = astronomical unit in cm */
/* MSUN = mass of the Sun in g */



/* Local */

/* ------------------------------------------------------------------------------ */

/* Change from longitude of perihelion to argument of perihelion */
    g = *p - *n;

/* Rotation factors */
    mco_sine__(i__, &si, &ci);
    mco_sine__(&g, &sg, &cg);
    mco_sine__(n, &sn, &cn);
    z1 = cg * cn;
    z2 = cg * sn;
    z3 = sg * cn;
    z4 = sg * sn;
    d11 = z1 - z4 * ci;
    d12 = z2 + z3 * ci;
    d13 = sg * si;
    d21 = -z3 - z2 * ci;
    d22 = -z4 + z1 * ci;
    d23 = cg * si;

/* Semi-major axis */
    a = *q / (1. - *e);

/* Ellipse */
    if (*e < 1.) {
	romes = sqrt(1. - *e * *e);
	temp = mco_kep__(e, l);
	mco_sine__(&temp, &se, &ce);
	z1 = a * (ce - *e);
	z2 = a * romes * se;
	temp = sqrt(*mu / a) / (1. - *e * ce);
	z3 = -se * temp;
	z4 = romes * ce * temp;
    } else {
/* Parabola */
	if (*e == 1.) {
	    ce = orbel_zget__(l);
	    z1 = *q * (1. - ce * ce);
	    z2 = *q * 2. * ce;
	    z4 = sqrt(*mu * 2. / *q) / (ce * ce + 1.);
	    z3 = -ce * z4;
	} else {
/* Hyperbola */
	    romes = sqrt(*e * *e - 1.);
	    temp = orbel_fhybrid__(e, l);
	    mco_sinh__(&temp, &se, &ce);
	    z1 = a * (ce - *e);
	    z2 = -a * romes * se;
	    temp = sqrt(*mu / abs(a)) / (*e * ce - 1.);
	    z3 = -se * temp;
	    z4 = romes * ce * temp;
	}
    }

    *x = d11 * z1 + d21 * z2;
    *y = d12 * z1 + d22 * z2;
    *z__ = d13 * z1 + d23 * z2;
    *u = d11 * z3 + d21 * z4;
    *v = d12 * z3 + d22 * z4;
    *w = d13 * z3 + d23 * z4;

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_el2x__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_KEP.FOR    (ErikSoft  7 July 1999) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Solves Kepler's equation for eccentricities less than one. */
/* Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330. */

/*  e = eccentricity */
/*  l = mean anomaly      (radians) */
/*  u = eccentric anomaly (   "   ) */

/* ------------------------------------------------------------------------------ */

doublereal mco_kep__(doublereal *e, doublereal *oldl)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double d_mod(doublereal *, doublereal *), sqrt(doublereal), log(
	    doublereal), exp(doublereal);

    /* Local variables */
    static doublereal l, p, q, x, f0, f1, f2, f3, p2, u1, u2, x2, z1, z2, z3, 
	    cc, pi, sn, ss;
    static logical big;
    static doublereal ome, dsn;
    static logical bigg, flag__;
    static doublereal sign, piby2, twopi;


/* Input/Outout */

/* Local */

/* ------------------------------------------------------------------------------ */

    pi = 3.141592653589793;
    twopi = pi * 2.;
    piby2 = pi * .5;

/* Reduce mean anomaly to lie in the range 0 < l < pi */
    if (*oldl >= 0.) {
	l = d_mod(oldl, &twopi);
    } else {
	l = d_mod(oldl, &twopi) + twopi;
    }
    sign = 1.;
    if (l > pi) {
	l = twopi - l;
	sign = -1.;
    }

    ome = 1. - *e;

    if (l >= .45 || *e < .55) {

/* Regions A,B or C in Nijenhuis */
/* ----------------------------- */

/* Rough starting value for eccentric anomaly */
	if (l < ome) {
	    u1 = ome;
	} else {
	    if (l > pi - 1. - *e) {
		u1 = (l + *e * pi) / (*e + 1.);
	    } else {
		u1 = l + *e;
	    }
	}

/* Improved value using Halley's method */
	flag__ = u1 > piby2;
	if (flag__) {
	    x = pi - u1;
	} else {
	    x = u1;
	}
	x2 = x * x;
	sn = x * (x2 * (x2 * .00761f - .16605f) + 1.);
	dsn = x2 * (x2 * .03805f - .49815f) + 1.;
	if (flag__) {
	    dsn = -dsn;
	}
	f2 = *e * sn;
	f0 = u1 - f2 - l;
	f1 = 1. - *e * dsn;
	u2 = u1 - f0 / (f1 - f0 * .5 * f2 / f1);
    } else {

/* Region D in Nijenhuis */
/* --------------------- */

/* Rough starting value for eccentric anomaly */
	z1 = *e * 4. + .5;
	p = ome / z1;
	q = l * .5 / z1;
	p2 = p * p;
	z2 = exp(log(sqrt(p2 * p + q * q) + q) / 1.5f);
	u1 = q * 2. / (z2 + p + p2 / z2);

/* Improved value using Newton's method */
	z2 = u1 * u1;
	z3 = z2 * z2;
	u2 = u1 - u1 * .075 * z3 / (ome + z1 * z2 + z3 * .375);
	u2 = l + *e * u2 * (3. - u2 * 4. * u2);
    }

/* Accurate value using 3rd-order version of Newton's method */
/* N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy! */

/* First get accurate values for u2 - sin(u2) and 1 - cos(u2) */
    bigg = u2 > piby2;
    if (bigg) {
	z3 = pi - u2;
    } else {
	z3 = u2;
    }

    big = z3 > piby2 * .5;
    if (big) {
	x = piby2 - z3;
    } else {
	x = z3;
    }

    x2 = x * x;
    ss = 1.;
    cc = 1.;

    ss = x * x2 / 6.f * (1.f - x2 / 20.f * (1.f - x2 / 42.f * (1.f - x2 / 
	    72.f * (1.f - x2 / 110.f * (1.f - x2 / 156.f * (1.f - x2 / 210.f *
	     (1.f - x2 / 272.f)))))));
    cc = x2 / 2.f * (1.f - x2 / 12.f * (1.f - x2 / 30.f * (1.f - x2 / 56.f * (
	    1.f - x2 / 90.f * (1.f - x2 / 132.f * (1.f - x2 / 182.f * (1.f - 
	    x2 / 240.f * (1.f - x2 / 306.f))))))));

    if (big) {
	z1 = cc + z3 - 1.;
	z2 = ss + z3 + 1. - piby2;
    } else {
	z1 = ss;
	z2 = cc;
    }

    if (bigg) {
	z1 = u2 * 2. + z1 - pi;
	z2 = 2. - z2;
    }

    f0 = l - u2 * ome - *e * z1;
    f1 = ome + *e * z2;
    f2 = *e * .5 * (u2 - z1);
    f3 = *e / 6. * (1. - z2);
    z1 = f0 / f1;
    z2 = f0 / (f2 * z1 + f1);
    ret_val = sign * (u2 + f0 / ((f3 * z1 + f2) * z2 + f1));

/* ------------------------------------------------------------------------------ */

    return ret_val;
} /* mco_kep__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_SINE.FOR    (ErikSoft  17 April 1997) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates sin and cos of an angle X (in radians). */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_sine__(doublereal *x, doublereal *sx, doublereal *cx)
{
    /* Builtin functions */
    double d_mod(doublereal *, doublereal *), cos(doublereal), sqrt(
	    doublereal);

    /* Local variables */
    static doublereal pi, twopi;



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    pi = 3.141592653589793;
    twopi = pi * 2.;

    if (*x > 0.) {
	*x = d_mod(x, &twopi);
    } else {
	*x = d_mod(x, &twopi) + twopi;
    }

    *cx = cos(*x);

    if (*x > pi) {
	*sx = -sqrt(1. - *cx * *cx);
    } else {
	*sx = sqrt(1. - *cx * *cx);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_sine__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_SINH.FOR    (ErikSoft  12 June 1998) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Calculates sinh and cosh of an angle X (in radians) */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_sinh__(doublereal *x, doublereal *sx, doublereal *cx)
{
    /* Builtin functions */
    double sinh(doublereal), sqrt(doublereal);



/* Input/Output */

/* ------------------------------------------------------------------------------ */

    *sx = sinh(*x);
    *cx = sqrt(*sx * *sx + 1.);

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_sinh__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_AEI.FOR    (ErikSoft   31 January 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Creates a filename and opens a file to store aei information for an object. */
/* The filename is based on the name of the object. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mio_aei__(char *id, char *extn, integer *unitnum, char *
	header, integer *lenhead, char *mem, integer *lmem, ftnlen id_len, 
	ftnlen extn_len, ftnlen header_len, ftnlen mem_len)
{
    /* Initialized data */

    static char bad[1*5] = "*" "/" "." ":" "&";

    /* System generated locals */
    integer i__1, i__2;
    olist o__1;
    inlist ioin__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_inqu(inlist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void), f_open(olist *);

    /* Local variables */
    static char filename[250];
    static integer j, k, lim[8]	/* was [2][4] */, nsub, itmp;
    static logical test;
    extern /* Subroutine */ int mio_spl__(integer *, char *, integer *, 
	    integer *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___177 = { 0, 6, 0, "(/,3a)", 0 };
    static cilist io___178 = { 0, 0, 0, "(/,30x,a8,//,a)", 0 };




/* Input/Output */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MERCURY.INC    (ErikSoft   4 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Parameters that you may want to alter at some point: */

/* NMAX  = maximum number of bodies */
/* CMAX  = maximum number of close-encounter minima monitored simultaneously */
/* NMESS = maximum number of messages in message.in */
/* HUGE  = an implausibly large number */
/* NFILES = maximum number of files that can be open at the same time */



/* ------------------------------------------------------------------------------ */

/* Constants: */

/* DR = conversion factor from degrees to radians */
/* K2 = Gaussian gravitational constant squared */
/* AU = astronomical unit in cm */
/* MSUN = mass of the Sun in g */



/* Local */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    --lmem;
    mem -= 80;

    /* Function Body */

/* Create a filename based on the object's name */
    mio_spl__(&c__8, id, &nsub, lim, (ftnlen)8);
/* Computing MIN */
    i__1 = 7, i__2 = lim[1] - lim[0];
    itmp = min(i__1,i__2);
    s_copy(filename, id, itmp + 1, itmp + 1);
    i__1 = itmp + 1;
    s_copy(filename + i__1, extn, itmp + 5 - i__1, (ftnlen)4);
    for (j = itmp + 6; j <= 250; ++j) {
	*(unsigned char *)&filename[j - 1] = ' ';
    }

/* Check for inappropriate characters in the filename */
    i__1 = itmp + 1;
    for (j = 1; j <= i__1; ++j) {
	for (k = 1; k <= 5; ++k) {
	    if (*(unsigned char *)&filename[j - 1] == *(unsigned char *)&bad[
		    k - 1]) {
		*(unsigned char *)&filename[j - 1] = '_';
	    }
	}
    }

/* If the file exists already, give a warning and don't overwrite it */
    ioin__1.inerr = 0;
    ioin__1.infilen = 250;
    ioin__1.infile = filename;
    ioin__1.inex = &test;
    ioin__1.inopen = 0;
    ioin__1.innum = 0;
    ioin__1.innamed = 0;
    ioin__1.inname = 0;
    ioin__1.inacc = 0;
    ioin__1.inseq = 0;
    ioin__1.indir = 0;
    ioin__1.infmt = 0;
    ioin__1.inform = 0;
    ioin__1.inunf = 0;
    ioin__1.inrecl = 0;
    ioin__1.innrec = 0;
    ioin__1.inblank = 0;
    f_inqu(&ioin__1);
    if (test) {
	s_wsfe(&io___177);
	do_fio(&c__1, mem + 9680, lmem[121]);
	do_fio(&c__1, mem + 6960, lmem[87]);
	do_fio(&c__1, filename, (ftnlen)80);
	e_wsfe();
	*unitnum = -1;
    } else {
	o__1.oerr = 0;
	o__1.ounit = *unitnum;
	o__1.ofnmlen = 250;
	o__1.ofnm = filename;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	io___178.ciunit = *unitnum;
	s_wsfe(&io___178);
	do_fio(&c__1, id, (ftnlen)8);
	do_fio(&c__1, header, (*lenhead));
	e_wsfe();
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mio_aei__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_C2FL.FOR    (ErikSoft   5 June 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* CHARACTER*8 ASCII string into a REAL*8 variable. */

/* N.B. X will lie in the range -1.e112 < X < 1.e112 */
/* === */

/* ------------------------------------------------------------------------------ */

doublereal mio_c2fl__(char *c__, ftnlen c_len)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal x;
    static integer ex;
    extern doublereal mio_c2re__(char *, doublereal *, doublereal *, integer *
	    , ftnlen);



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    x = mio_c2re__(c__, &c_b104, &c_b118, &c__7, (ftnlen)8);
    x = x * 2. - 1.;
    ex = (*(unsigned char *)&c__[7] + 256) % 256 - 144;
    d__1 = (doublereal) ex;
    ret_val = x * pow_dd(&c_b198, &d__1);

/* ------------------------------------------------------------------------------ */

    return ret_val;
} /* mio_c2fl__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_C2RE.FOR    (ErikSoft   5 June 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Converts an ASCII string into a REAL*8 variable X, where XMIN <= X < XMAX, */
/* using the new format compression: */

/* X is assumed to be made up of NCHAR base-224 digits, each one represented */
/* by a character in the ASCII string. Each digit is given by the ASCII */
/* number of the character minus 32. */
/* The first 32 ASCII characters (CTRL characters) are avoided, because they */
/* cause problems when using some operating systems. */

/* ------------------------------------------------------------------------------ */

doublereal mio_c2re__(char *c__, doublereal *xmin, doublereal *xmax, integer *
	nchar, ftnlen c_len)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer j;
    static doublereal y;



/* Input/output */

/* Local */

/* ------------------------------------------------------------------------------ */

    y = 0.;
    for (j = *nchar; j >= 1; --j) {
	y = (y + (doublereal) ((*(unsigned char *)&c__[j - 1] + 256) % 256 - 
		32)) / 224.;
    }

    ret_val = *xmin + y * (*xmax - *xmin);

/* ------------------------------------------------------------------------------ */

    return ret_val;
} /* mio_c2re__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_ERR.FOR    (ErikSoft  6 December 1999) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Writes out an error message and terminates Mercury. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mio_err__(integer *unit, char *s1, integer *ls1, char *
	s2, integer *ls2, char *s3, integer *ls3, char *s4, integer *ls4, 
	ftnlen s1_len, ftnlen s2_len, ftnlen s3_len, ftnlen s4_len)
{
    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___183 = { 0, 6, 0, "(a)", 0 };
    static cilist io___184 = { 0, 0, 0, "(/,3a,/,2a)", 0 };




/* Input/Output */

/* ------------------------------------------------------------------------------ */

    s_wsfe(&io___183);
    do_fio(&c__1, " ERROR: Programme terminated.", (ftnlen)29);
    e_wsfe();
    io___184.ciunit = *unit;
    s_wsfe(&io___184);
    do_fio(&c__1, s1, (*ls1));
    do_fio(&c__1, s2, (*ls2));
    do_fio(&c__1, s3, (*ls3));
    do_fio(&c__1, " ", (ftnlen)1);
    do_fio(&c__1, s4, (*ls4));
    e_wsfe();
    s_stop("", (ftnlen)0);

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mio_err__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_H2B.FOR    (ErikSoft   2 November 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Converts coordinates with respect to the central body to barycentric */
/* coordinates. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_h2b__(doublereal *jcen, integer *nbod, integer *nbig,
	 doublereal *h__, doublereal *m, doublereal *xh, doublereal *vh, 
	doublereal *x, doublereal *v)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal temp, mtot;



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    --jcen;
    v -= 4;
    x -= 4;
    vh -= 4;
    xh -= 4;
    --m;

    /* Function Body */
    mtot = 0.;
    x[4] = 0.;
    x[5] = 0.;
    x[6] = 0.;
    v[4] = 0.;
    v[5] = 0.;
    v[6] = 0.;

/* Calculate coordinates and velocities of the central body */
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	mtot += m[j];
	x[4] += m[j] * xh[j * 3 + 1];
	x[5] += m[j] * xh[j * 3 + 2];
	x[6] += m[j] * xh[j * 3 + 3];
	v[4] += m[j] * vh[j * 3 + 1];
	v[5] += m[j] * vh[j * 3 + 2];
	v[6] += m[j] * vh[j * 3 + 3];
    }

    temp = -1. / (mtot + m[1]);
    x[4] = temp * x[4];
    x[5] = temp * x[5];
    x[6] = temp * x[6];
    v[4] = temp * v[4];
    v[5] = temp * v[5];
    v[6] = temp * v[6];

/* Calculate the barycentric coordinates and velocities */
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	x[j * 3 + 1] = xh[j * 3 + 1] + x[4];
	x[j * 3 + 2] = xh[j * 3 + 2] + x[5];
	x[j * 3 + 3] = xh[j * 3 + 3] + x[6];
	v[j * 3 + 1] = vh[j * 3 + 1] + v[4];
	v[j * 3 + 2] = vh[j * 3 + 2] + v[5];
	v[j * 3 + 3] = vh[j * 3 + 3] + v[6];
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_h2b__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_H2CB.FOR    (ErikSoft   2 November 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Convert coordinates with respect to the central body to close-binary */
/* coordinates. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_h2cb__(doublereal *jcen, integer *nbod, integer *
	nbig, doublereal *h__, doublereal *m, doublereal *xh, doublereal *vh, 
	doublereal *x, doublereal *v)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal mbin, temp, msum, mbin_1__, mvsum[3], mtot_1__;



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    --jcen;
    v -= 4;
    x -= 4;
    vh -= 4;
    xh -= 4;
    --m;

    /* Function Body */
    msum = 0.;
    mvsum[0] = 0.;
    mvsum[1] = 0.;
    mvsum[2] = 0.;
    mbin = m[1] + m[2];
    mbin_1__ = 1. / mbin;

    x[7] = xh[7];
    x[8] = xh[8];
    x[9] = xh[9];
    temp = m[1] * mbin_1__;
    v[7] = temp * vh[7];
    v[8] = temp * vh[8];
    v[9] = temp * vh[9];

    i__1 = *nbod;
    for (j = 3; j <= i__1; ++j) {
	msum += m[j];
	mvsum[0] += m[j] * vh[j * 3 + 1];
	mvsum[1] += m[j] * vh[j * 3 + 2];
	mvsum[2] += m[j] * vh[j * 3 + 3];
    }
    mtot_1__ = 1. / (msum + mbin);
    mvsum[0] = mtot_1__ * (mvsum[0] + m[2] * vh[7]);
    mvsum[1] = mtot_1__ * (mvsum[1] + m[2] * vh[8]);
    mvsum[2] = mtot_1__ * (mvsum[2] + m[2] * vh[9]);

    temp = m[2] * mbin_1__;
    i__1 = *nbod;
    for (j = 3; j <= i__1; ++j) {
	x[j * 3 + 1] = xh[j * 3 + 1] - temp * xh[7];
	x[j * 3 + 2] = xh[j * 3 + 2] - temp * xh[8];
	x[j * 3 + 3] = xh[j * 3 + 3] - temp * xh[9];
	v[j * 3 + 1] = vh[j * 3 + 1] - mvsum[0];
	v[j * 3 + 2] = vh[j * 3 + 2] - mvsum[1];
	v[j * 3 + 3] = vh[j * 3 + 3] - mvsum[2];
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_h2cb__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_H2J.FOR    (ErikSoft   2 November 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Converts coordinates with respect to the central body to Jacobi coordinates. */
/* Note that the Jacobi coordinates of all small bodies are assumed to be the */
/* same as their coordinates with respect to the central body. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_h2j__(doublereal *jcen, integer *nbod, integer *nbig,
	 doublereal *h__, doublereal *m, doublereal *xh, doublereal *vh, 
	doublereal *x, doublereal *v)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal mu, mv, mw, mx, my, mz, temp, mtot;



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------c */
    /* Parameter adjustments */
    --jcen;
    v -= 4;
    x -= 4;
    vh -= 4;
    xh -= 4;
    --m;

    /* Function Body */
    mtot = m[2];
    x[7] = xh[7];
    x[8] = xh[8];
    x[9] = xh[9];
    v[7] = vh[7];
    v[8] = vh[8];
    v[9] = vh[9];
    mx = m[2] * xh[7];
    my = m[2] * xh[8];
    mz = m[2] * xh[9];
    mu = m[2] * vh[7];
    mv = m[2] * vh[8];
    mw = m[2] * vh[9];

    i__1 = *nbig - 1;
    for (j = 3; j <= i__1; ++j) {
	temp = 1. / (mtot + m[1]);
	mtot += m[j];
	x[j * 3 + 1] = xh[j * 3 + 1] - temp * mx;
	x[j * 3 + 2] = xh[j * 3 + 2] - temp * my;
	x[j * 3 + 3] = xh[j * 3 + 3] - temp * mz;
	v[j * 3 + 1] = vh[j * 3 + 1] - temp * mu;
	v[j * 3 + 2] = vh[j * 3 + 2] - temp * mv;
	v[j * 3 + 3] = vh[j * 3 + 3] - temp * mw;
	mx += m[j] * xh[j * 3 + 1];
	my += m[j] * xh[j * 3 + 2];
	mz += m[j] * xh[j * 3 + 3];
	mu += m[j] * vh[j * 3 + 1];
	mv += m[j] * vh[j * 3 + 2];
	mw += m[j] * vh[j * 3 + 3];
    }

    if (*nbig > 2) {
	temp = 1. / (mtot + m[1]);
	x[*nbig * 3 + 1] = xh[*nbig * 3 + 1] - temp * mx;
	x[*nbig * 3 + 2] = xh[*nbig * 3 + 2] - temp * my;
	x[*nbig * 3 + 3] = xh[*nbig * 3 + 3] - temp * mz;
	v[*nbig * 3 + 1] = vh[*nbig * 3 + 1] - temp * mu;
	v[*nbig * 3 + 2] = vh[*nbig * 3 + 2] - temp * mv;
	v[*nbig * 3 + 3] = vh[*nbig * 3 + 3] - temp * mw;
    }

    i__1 = *nbod;
    for (j = *nbig + 1; j <= i__1; ++j) {
	x[j * 3 + 1] = xh[j * 3 + 1];
	x[j * 3 + 2] = xh[j * 3 + 2];
	x[j * 3 + 3] = xh[j * 3 + 3];
	v[j * 3 + 1] = vh[j * 3 + 1];
	v[j * 3 + 2] = vh[j * 3 + 2];
	v[j * 3 + 3] = vh[j * 3 + 3];
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_h2j__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_IDEN.FOR    (ErikSoft   2 November 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Makes a new copy of a set of coordinates. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_iden__(doublereal *jcen, integer *nbod, integer *
	nbig, doublereal *h__, doublereal *m, doublereal *xh, doublereal *vh, 
	doublereal *x, doublereal *v)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    --jcen;
    v -= 4;
    x -= 4;
    vh -= 4;
    xh -= 4;
    --m;

    /* Function Body */
    i__1 = *nbod;
    for (j = 1; j <= i__1; ++j) {
	x[j * 3 + 1] = xh[j * 3 + 1];
	x[j * 3 + 2] = xh[j * 3 + 2];
	x[j * 3 + 3] = xh[j * 3 + 3];
	v[j * 3 + 1] = vh[j * 3 + 1];
	v[j * 3 + 2] = vh[j * 3 + 2];
	v[j * 3 + 3] = vh[j * 3 + 3];
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_iden__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_X2EL.FOR    (ErikSoft  20 February 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates Keplerian orbital elements given relative coordinates and */
/* velocities, and GM = G times the sum of the masses. */

/* The elements are: q = perihelion distance */
/*                   e = eccentricity */
/*                   i = inclination */
/*                   p = longitude of perihelion (NOT argument of perihelion!!) */
/*                   n = longitude of ascending node */
/*                   l = mean anomaly (or mean longitude if e < 1.e-8) */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_x2el__(doublereal *gm, doublereal *x, doublereal *y, 
	doublereal *z__, doublereal *u, doublereal *v, doublereal *w, 
	doublereal *q, doublereal *e, doublereal *i__, doublereal *p, 
	doublereal *n, doublereal *l)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal), atan2(doublereal, doublereal), 
	    d_sign(doublereal *, doublereal *), sin(doublereal), log(
	    doublereal), sinh(doublereal), d_mod(doublereal *, doublereal *);

    /* Local variables */
    static doublereal f, h__, r__, s, h2, v2, ce, cf, ci, hx, hy, hz, to, rv, 
	    tmp2, bige, temp, true__;



/* Input/Output */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MERCURY.INC    (ErikSoft   4 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Parameters that you may want to alter at some point: */

/* NMAX  = maximum number of bodies */
/* CMAX  = maximum number of close-encounter minima monitored simultaneously */
/* NMESS = maximum number of messages in message.in */
/* HUGE  = an implausibly large number */
/* NFILES = maximum number of files that can be open at the same time */



/* ------------------------------------------------------------------------------ */

/* Constants: */

/* DR = conversion factor from degrees to radians */
/* K2 = Gaussian gravitational constant squared */
/* AU = astronomical unit in cm */
/* MSUN = mass of the Sun in g */



/* Local */

/* ------------------------------------------------------------------------------ */

    hx = *y * *w - *z__ * *v;
    hy = *z__ * *u - *x * *w;
    hz = *x * *v - *y * *u;
    h2 = hx * hx + hy * hy + hz * hz;
    v2 = *u * *u + *v * *v + *w * *w;
    rv = *x * *u + *y * *v + *z__ * *w;
    r__ = sqrt(*x * *x + *y * *y + *z__ * *z__);
    h__ = sqrt(h2);
    s = h2 / *gm;

/* Inclination and node */
    ci = hz / h__;
    if (abs(ci) < 1.) {
	*i__ = acos(ci);
	*n = atan2(hx, -hy);
	if (*n < 0.) {
	    *n += 6.2831853071795862;
	}
    } else {
	if (ci > 0.) {
	    *i__ = 0.;
	}
	if (ci < 0.) {
	    *i__ = 3.141592653589793;
	}
	*n = 0.;
    }

/* Eccentricity and perihelion distance */
    temp = s * (v2 / *gm - 2. / r__) + 1.;
    if (temp <= 0.) {
	*e = 0.;
    } else {
	*e = sqrt(temp);
    }
    *q = s / (*e + 1.);

/* True longitude */
    if (hy != 0.) {
	to = -hx / hy;
	temp = (1. - ci) * to;
	tmp2 = to * to;
	true__ = atan2(*y * (tmp2 * ci + 1.) - *x * temp, *x * (tmp2 + ci) - *
		y * temp);
    } else {
	true__ = atan2(*y * ci, *x);
    }
    if (ci < 0.) {
	true__ += 3.141592653589793;
    }

    if (*e < 3e-8) {
	*p = 0.;
	*l = true__;
    } else {
	ce = (v2 * r__ - *gm) / (*e * *gm);

/* Mean anomaly for ellipse */
	if (*e < 1.) {
	    if (abs(ce) > 1.) {
		ce = d_sign(&c_b118, &ce);
	    }
	    bige = acos(ce);
	    if (rv < 0.) {
		bige = 6.2831853071795862 - bige;
	    }
	    *l = bige - *e * sin(bige);
	} else {

/* Mean anomaly for hyperbola */
	    if (ce < 1.) {
		ce = 1.;
	    }
	    bige = log(ce + sqrt(ce * ce - 1.));
	    if (rv < 0.) {
		bige = -bige;
	    }
	    *l = *e * sinh(bige) - bige;
	}

/* Longitude of perihelion */
	cf = (s - r__) / (*e * r__);
	if (abs(cf) > 1.) {
	    cf = d_sign(&c_b118, &cf);
	}
	f = acos(cf);
	if (rv < 0.) {
	    f = 6.2831853071795862 - f;
	}
	*p = true__ - f;
	d__1 = *p + 6.2831853071795862 + 6.2831853071795862;
	*p = d_mod(&d__1, &c_b141);
    }

    if (*l < 0. && *e < 1.) {
	*l += 6.2831853071795862;
    }
    if (*l > 6.2831853071795862 && *e < 1.) {
	*l = d_mod(l, &c_b141);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_x2el__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_JD_Y.FOR    (ErikSoft  2 June 1998) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Converts from Julian day number to Julian/Gregorian Calendar dates, assuming */
/* the dates are those used by the English calendar. */

/* Algorithm taken from `Practical Astronomy with your calculator' (1988) */
/* by Peter Duffett-Smith, 3rd edition, C.U.P. */

/* Algorithm for negative Julian day numbers (Julian calendar assumed) by */
/* J. E. Chambers. */

/* N.B. The output date is with respect to the Julian Calendar on or before */
/* ===  4th October 1582, and with respect to the Gregorian Calendar on or */
/*      after 15th October 1582. */


/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mio_jd_y__(doublereal *jd0, integer *year, integer *
	month, doublereal *day)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_int(doublereal *), d_sign(doublereal *, doublereal *), d_mod(
	    doublereal *, doublereal *);

    /* Local variables */
    static integer a, b, c__, d__, e;
    static doublereal f;
    static integer g, i__;
    static doublereal x, y, z__, jd, temp;



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    if (*jd0 <= 0.) {
	goto L50;
    }

    jd = *jd0 + .5;
    d__2 = abs(jd);
    d__1 = d_int(&d__2);
    i__ = (integer) d_sign(&d__1, &jd);
    f = jd - i__ * 1.;

/* If on or after 15th October 1582 */
    if (i__ > 2299160) {
	temp = (i__ * 1. - 1867216.25) / 36524.25;
	d__2 = abs(temp);
	d__1 = d_int(&d__2);
	a = (integer) d_sign(&d__1, &temp);
	temp = a * .25;
	d__2 = abs(temp);
	d__1 = d_int(&d__2);
	b = (integer) (i__ + 1 + a - d_sign(&d__1, &temp));
    } else {
	b = i__;
    }

    c__ = b + 1524;
    temp = (c__ * 1. - 122.1) / 365.25;
    d__2 = abs(temp);
    d__1 = d_int(&d__2);
    d__ = (integer) d_sign(&d__1, &temp);
    temp = d__ * 365.25;
    d__2 = abs(temp);
    d__1 = d_int(&d__2);
    e = (integer) d_sign(&d__1, &temp);
    temp = (c__ - e) / 30.6001;
    d__2 = abs(temp);
    d__1 = d_int(&d__2);
    g = (integer) d_sign(&d__1, &temp);

    temp = g * 30.6001;
    d__2 = abs(temp);
    d__1 = d_int(&d__2);
    *day = (c__ - e) * 1. + f - d_sign(&d__1, &temp) * 1.;

    if (g <= 13) {
	*month = g - 1;
    }
    if (g > 13) {
	*month = g - 13;
    }

    if (*month > 2) {
	*year = d__ - 4716;
    }
    if (*month <= 2) {
	*year = d__ - 4715;
    }

    if (*day > 32.) {
	*day += -32;
	++(*month);
    }

    if (*month > 12) {
	*month += -12;
	++(*year);
    }
    return 0;

L50:

/* Algorithm for negative Julian day numbers (Duffett-Smith won't work) */
    x = *jd0 - 2232101.5f;
    f = x - d_int(&x);
    if (f < 0.) {
	f += 1.;
    }
    d__1 = d_mod(&x, &c_b245) + 1461.;
    y = d_int(&d__1);
    d__1 = d_mod(&y, &c_b246);
    z__ = d_int(&d__1);
    *month = (integer) ((z__ + .5) / 30.61);
    d__1 = z__ + 1.5 - (doublereal) (*month) * 30.61;
    *day = d_int(&d__1) + f;
    *month = (*month + 2) % 12 + 1;

    *year = (integer) (x / 365.25) + 1399;
    if (x < 0.) {
	--(*year);
    }
    if (*month < 3) {
	++(*year);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mio_jd_y__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_SPL.FOR    (ErikSoft  14 November 1999) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Given a character string STRING, of length LEN bytes, the routine finds */
/* the beginnings and ends of NSUB substrings present in the original, and */
/* delimited by spaces. The positions of the extremes of each substring are */
/* returned in the array DELIMIT. */
/* Substrings are those which are separated by spaces or the = symbol. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mio_spl__(integer *len, char *string, integer *nsub, 
	integer *delimit, ftnlen string_len)
{
    static char c__[1];
    static integer j, k;



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    --string;
    delimit -= 3;

    /* Function Body */
    *nsub = 0;
    j = 0;
    *(unsigned char *)c__ = ' ';
    delimit[3] = -1;

/* Find the start of string */
L10:
    ++j;
    if (j > *len) {
	goto L99;
    }
    *(unsigned char *)c__ = *(unsigned char *)&string[j];
    if (*(unsigned char *)c__ == ' ' || *(unsigned char *)c__ == '=') {
	goto L10;
    }

/* Find the end of string */
    k = j;
L20:
    ++k;
    if (k > *len) {
	goto L30;
    }
    *(unsigned char *)c__ = *(unsigned char *)&string[k];
    if (*(unsigned char *)c__ != ' ' && *(unsigned char *)c__ != '=') {
	goto L20;
    }

/* Store details for this string */
L30:
    ++(*nsub);
    delimit[(*nsub << 1) + 1] = j;
    delimit[(*nsub << 1) + 2] = k - 1;

    if (k < *len) {
	j = k;
	goto L10;
    }

L99:

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mio_spl__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MXX_SORT.FOR    (ErikSoft 24 May 1997) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Sorts an array X, of size N, using Shell's method. Also returns an array */
/* INDEX that gives the original index of every item in the sorted array X. */

/* N.B. The maximum array size is 29523. */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mxx_sort__(integer *n, doublereal *x, integer *index)
{
    /* Initialized data */

    static integer incarr[9] = { 1,4,13,40,121,364,1093,3280,9841 };

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal y;
    static integer iy, inc;



/* Input/Output */

/* Local */
    /* Parameter adjustments */
    --index;
    --x;

    /* Function Body */

/* ------------------------------------------------------------------------------ */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	index[i__] = i__;
    }

    m = 0;
L10:
    ++m;
    if (incarr[m - 1] < *n) {
	goto L10;
    }
    --m;

    for (i__ = m; i__ >= 1; --i__) {
	inc = incarr[i__ - 1];
	i__1 = inc;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *n - j;
	    i__3 = inc;
	    for (k = inc; i__3 < 0 ? k >= i__2 : k <= i__2; k += i__3) {
		y = x[j + k];
		iy = index[j + k];
		i__4 = j;
		i__5 = -inc;
		for (l = j + k - inc; i__5 < 0 ? l >= i__4 : l <= i__4; l += 
			i__5) {
		    if (x[l] <= y) {
			goto L20;
		    }
		    x[l + inc] = x[l];
		    index[l + inc] = index[l];
		}
L20:
		x[l + inc] = y;
		index[l + inc] = iy;
	    }
	}
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mxx_sort__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      M_SFUNC.FOR     (ErikSoft  14 November 1998) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Calculates Z = [ (3 + S^2)arctan(S) - 3S ] / S^3 and its derivative DZ, */
/* for S > 0. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int m_sfunc__(doublereal *s, doublereal *z__, doublereal *dz)
{
    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static doublereal a, s2, s4, s6, s8;



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    s2 = *s * *s;

    if (*s > .01) {
	a = atan(*s);
	*z__ = ((s2 + 3.) * a - *s * 3.) / (*s * s2);
	*dz = (*s * 2. * a - 3. + (s2 + 3.) / (s2 + 1.)) / (*s * s2) - *z__ * 
		3. / *s;
    } else {
	s4 = s2 * s2;
	s6 = s2 * s4;
	s8 = s4 * s4;
	*z__ = s8 * -.1616161616161616 + s6 * .1904761904761905 - s4 * 
		.2285714285714286 + s2 * .2666666666666667;
	*dz = *s * (s6 * -1.292929292929293 + s4 * 1.142857142857143 - s2 * 
		.914285714285714 + .533333333333333);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* m_sfunc__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      M_FORMAT.FOR    (ErikSoft   31 January 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Makes an output format list and file header for the orbital-element files */
/* created by M_ELEM3.FOR */
/* Also identifies which orbital elements will be output for each object. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int m_format__(char *string, integer *timestyle, integer *
	nel, integer *iel, char *fout, char *header, integer *lenhead, ftnlen 
	string_len, ftnlen fout_len, ftnlen header_len)
{
    /* Initialized data */

    static char elcode[1*22] = "a" "e" "i" "g" "n" "l" "p" "q" "b" "x" "y" 
	    "z" "u" "v" "w" "r" "f" "m" "o" "s" "d" "c";
    static char elhead[4*22] = "  a " "  e " "  i " "peri" "node" "  M " 
	    "long" "  q " "  Q " "  x " "  y " "  z " " vx " " vy " " vz " 
	    "  r " "  f " "mass" "oblq" "spin" "dens" "comp";

    /* System generated locals */
    integer i__1, i__2, i__3;
    icilist ici__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void), s_wsfi(
	    icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void);

    /* Local variables */
    static integer formflag, i__, j, f1, f2, lim[40]	/* was [2][20] */, 
	    pos, nsub, itmp;
    extern /* Subroutine */ int mio_spl__(integer *, char *, integer *, 
	    integer *, ftnlen);
    static integer lenfout;



/* Input/Output */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MERCURY.INC    (ErikSoft   4 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Parameters that you may want to alter at some point: */

/* NMAX  = maximum number of bodies */
/* CMAX  = maximum number of close-encounter minima monitored simultaneously */
/* NMESS = maximum number of messages in message.in */
/* HUGE  = an implausibly large number */
/* NFILES = maximum number of files that can be open at the same time */



/* ------------------------------------------------------------------------------ */

/* Constants: */

/* DR = conversion factor from degrees to radians */
/* K2 = Gaussian gravitational constant squared */
/* AU = astronomical unit in cm */
/* MSUN = mass of the Sun in g */



/* Local */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    --iel;

    /* Function Body */

/* Initialize header to a blank string */
    for (i__ = 1; i__ <= 250; ++i__) {
	*(unsigned char *)&header[i__ - 1] = ' ';
    }

/* Create part of the format list and header for the required time style */
    if (*timestyle == 0 || *timestyle == 2) {
	s_copy(fout, "(1x,f18.5", (ftnlen)9, (ftnlen)9);
	lenfout = 9;
	s_copy(header, "    Time (days)    ", (ftnlen)19, (ftnlen)19);
	*lenhead = 19;
    } else if (*timestyle == 1) {
	s_copy(fout, "(1x,i10,1x,i2,1x,f8.5", (ftnlen)21, (ftnlen)21);
	lenfout = 21;
	s_copy(header, "    Year/Month/Day     ", (ftnlen)23, (ftnlen)23);
	*lenhead = 23;
    } else if (*timestyle == 3) {
	s_copy(fout, "(1x,f18.7", (ftnlen)9, (ftnlen)9);
	lenfout = 9;
	s_copy(header, "    Time (years)   ", (ftnlen)19, (ftnlen)19);
	*lenhead = 19;
    }

/* Identify the required elements */
    mio_spl__(&c__250, string, &nsub, lim, (ftnlen)250);
    i__1 = nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
	for (j = 1; j <= 22; ++j) {
	    i__2 = lim[(i__ << 1) - 2] - 1;
	    if (s_cmp(string + i__2, elcode + (j - 1), lim[(i__ << 1) - 2] - 
		    i__2, (ftnlen)1) == 0) {
		iel[i__] = j;
	    }
	}
    }
    *nel = nsub;

/* For each element, see whether normal or exponential notation is required */
    i__1 = nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
	formflag = 0;
	i__2 = lim[(i__ << 1) - 1];
	for (j = lim[(i__ << 1) - 2] + 1; j <= i__2; ++j) {
	    if (formflag == 0) {
		pos = j;
	    }
	    if (*(unsigned char *)&string[j - 1] == '.') {
		formflag = 1;
	    }
	    if (*(unsigned char *)&string[j - 1] == 'e') {
		formflag = 2;
	    }
	}

/* Create the rest of the format list and header */
	if (formflag == 1) {
	    i__2 = lim[(i__ << 1) - 2];
	    ici__1.icierr = 0;
	    ici__1.iciend = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = pos - 1 - i__2;
	    ici__1.iciunit = string + i__2;
	    ici__1.icifmt = 0;
	    s_rsli(&ici__1);
	    do_lio(&c__3, &c__1, (char *)&f1, (ftnlen)sizeof(integer));
	    e_rsli();
	    i__2 = pos;
	    ici__1.icierr = 0;
	    ici__1.iciend = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = lim[(i__ << 1) - 1] - i__2;
	    ici__1.iciunit = string + i__2;
	    ici__1.icifmt = 0;
	    s_rsli(&ici__1);
	    do_lio(&c__3, &c__1, (char *)&f2, (ftnlen)sizeof(integer));
	    e_rsli();
	    i__2 = lenfout;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = lenfout + 10 - i__2;
	    ici__1.iciunit = fout + i__2;
	    ici__1.icifmt = "(a10)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, ",1x,f  .  ", (ftnlen)10);
	    e_wsfi();
	    i__2 = lenfout + 5;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = lenfout + 7 - i__2;
	    ici__1.iciunit = fout + i__2;
	    ici__1.icifmt = "(i2)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, (char *)&f1, (ftnlen)sizeof(integer));
	    e_wsfi();
	    i__2 = lenfout + 8;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = lenfout + 10 - i__2;
	    ici__1.iciunit = fout + i__2;
	    ici__1.icifmt = "(i2)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, (char *)&f2, (ftnlen)sizeof(integer));
	    e_wsfi();
	    lenfout += 10;
	} else if (formflag == 2) {
	    i__2 = lim[(i__ << 1) - 2];
	    ici__1.icierr = 0;
	    ici__1.iciend = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = pos - 1 - i__2;
	    ici__1.iciunit = string + i__2;
	    ici__1.icifmt = 0;
	    s_rsli(&ici__1);
	    do_lio(&c__3, &c__1, (char *)&f1, (ftnlen)sizeof(integer));
	    e_rsli();
	    i__2 = lenfout;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = lenfout + 16 - i__2;
	    ici__1.iciunit = fout + i__2;
	    ici__1.icifmt = "(a16)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, ",1x,1p,e  .  ,0p", (ftnlen)16);
	    e_wsfi();
	    i__2 = lenfout + 8;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = lenfout + 10 - i__2;
	    ici__1.iciunit = fout + i__2;
	    ici__1.icifmt = "(i2)";
	    s_wsfi(&ici__1);
	    do_fio(&c__1, (char *)&f1, (ftnlen)sizeof(integer));
	    e_wsfi();
	    i__2 = lenfout + 11;
	    ici__1.icierr = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = lenfout + 13 - i__2;
	    ici__1.iciunit = fout + i__2;
	    ici__1.icifmt = "(i2)";
	    s_wsfi(&ici__1);
	    i__3 = f1 - 7;
	    do_fio(&c__1, (char *)&i__3, (ftnlen)sizeof(integer));
	    e_wsfi();
	    lenfout += 16;
	}
	itmp = (f1 - 4) / 2;
	i__2 = *lenhead + itmp + 1;
	s_copy(header + i__2, elhead + (iel[i__] - 1 << 2), *lenhead + itmp + 
		5 - i__2, (ftnlen)4);
	*lenhead = *lenhead + f1 + 1;
    }

    ++lenfout;
    *(unsigned char *)&fout[lenfout - 1] = ')';

/* ------------------------------------------------------------------------------ */

    return 0;
} /* m_format__ */


/* ********************************************************************** */
/*                    ORBEL_FHYBRID.F */
/* ********************************************************************** */
/*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach. */

/*             Input: */
/*                           e ==> eccentricity anomaly. (real scalar) */
/*                           n ==> hyperbola mean anomaly. (real scalar) */
/*             Returns: */
/*               orbel_fhybrid ==>  eccentric anomaly. (real scalar) */

/*     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON */
/* 	         For larger N, uses FGET */
/*     REMARKS: */
/*     AUTHOR: M. Duncan */
/*     DATE WRITTEN: May 26,1992. */
/*     REVISIONS: */
/*     REVISIONS: 2/26/93 hfl */
/* ********************************************************************** */
doublereal orbel_fhybrid__(doublereal *e, doublereal *n)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    extern doublereal orbel_fget__(doublereal *, doublereal *), orbel_flon__(
	    doublereal *, doublereal *);
    static doublereal abn;

/* ...  Inputs Only: */
/* ************************************************************************* */
/*                        SWIFT.INC */
/* ************************************************************************* */
/* Include file for SWIFT */

/* Author:  Hal Levison */
/* Date:    2/2/93 */
/* Last revision: 3/7/93 */
/* ...   Maximum array size */
/* you got it baby */
/* max number of planets, including th */
/* ...   Size of the test particle status flag */
/* max number of test particles */
/* ...   convergence criteria for danby */
/* ...    loop limits in the Laguerre attempts */
/* ...    A small number */
/* ...    trig stuff */
/* ------------------------------------------------------------------------- */
/* ...  Internals: */
/* ---- */
/* ...  Executable code */
    abn = *n;
    if (*n < 0.) {
	abn = -abn;
    }
    if (abn < *e * .636 - .6) {
	ret_val = orbel_flon__(e, n);
    } else {
	ret_val = orbel_fget__(e, n);
    }
    return ret_val;
} /* orbel_fhybrid__ */

/* ------------------------------------------------------------------- */

/* ********************************************************************** */
/*                    ORBEL_FGET.F */
/* ********************************************************************** */
/*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach. */

/*             Input: */
/*                           e ==> eccentricity anomaly. (real scalar) */
/*                        capn ==> hyperbola mean anomaly. (real scalar) */
/*             Returns: */
/*                  orbel_fget ==>  eccentric anomaly. (real scalar) */

/*     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of */
/*           Cel. Mech. ".  Quartic convergence from Danby's book. */
/*     REMARKS: */
/*     AUTHOR: M. Duncan */
/*     DATE WRITTEN: May 11, 1992. */
/*     REVISIONS: 2/26/93 hfl */
/*     Modified by JEC */
/* ********************************************************************** */
/* orbel_fhybrid */
doublereal orbel_fget__(doublereal *e, doublereal *capn)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    extern /* Subroutine */ int mco_sinh__(doublereal *, doublereal *, 
	    doublereal *);
    static doublereal f;
    static integer i__;
    static doublereal x, fp, dx, ech, esh, chx, fpp, tmp, shx, fppp;

    /* Fortran I/O blocks */
    static cilist io___278 = { 0, 6, 0, 0, 0 };


/* ...  Inputs Only: */
/* ************************************************************************* */
/*                        SWIFT.INC */
/* ************************************************************************* */
/* Include file for SWIFT */

/* Author:  Hal Levison */
/* Date:    2/2/93 */
/* Last revision: 3/7/93 */
/* ...   Maximum array size */
/* you got it baby */
/* max number of planets, including th */
/* ...   Size of the test particle status flag */
/* max number of test particles */
/* ...   convergence criteria for danby */
/* ...    loop limits in the Laguerre attempts */
/* ...    A small number */
/* ...    trig stuff */
/* ------------------------------------------------------------------------- */
/* ...  Internals: */
/* ---- */
/* ...  Executable code */
/* Function to solve "Kepler's eqn" for F (here called */
/* x) for given e and CAPN. */
/*  begin with a guess proposed by Danby */
    if (*capn < 0.) {
	tmp = *capn * -2. / *e + 1.8;
	x = -log(tmp);
    } else {
	tmp = *capn * 2. / *e + 1.8;
	x = log(tmp);
    }
    ret_val = x;
    for (i__ = 1; i__ <= 10; ++i__) {
	mco_sinh__(&x, &shx, &chx);
	esh = *e * shx;
	ech = *e * chx;
	f = esh - x - *capn;
/* 	  write(6,*) 'i,x,f : ',i,x,f */
	fp = ech - 1.;
	fpp = esh;
	fppp = ech;
	dx = -f / fp;
	dx = -f / (fp + dx * fpp / 2.);
	dx = -f / (fp + dx * fpp / 2. + dx * dx * fppp / 6.);
	ret_val = x + dx;
/*   If we have converged here there's no point in going on */
	if (abs(dx) <= 4e-15) {
	    return ret_val;
	}
	x = ret_val;
    }
    s_wsle(&io___278);
    do_lio(&c__9, &c__1, "FGET : RETURNING WITHOUT COMPLETE CONVERGENCE", (
	    ftnlen)45);
    e_wsle();
    return ret_val;
} /* orbel_fget__ */

/* ------------------------------------------------------------------ */

/* ********************************************************************** */
/*                    ORBEL_FLON.F */
/* ********************************************************************** */
/*     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach. */

/*             Input: */
/*                           e ==> eccentricity anomaly. (real scalar) */
/*                        capn ==> hyperbola mean anomaly. (real scalar) */
/*             Returns: */
/*                  orbel_flon ==>  eccentric anomaly. (real scalar) */

/*     ALGORITHM: Uses power series for N in terms of F and Newton,s method */
/*     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6) */
/*     AUTHOR: M. Duncan */
/*     DATE WRITTEN: May 26, 1992. */
/*     REVISIONS: */
/* ********************************************************************** */
/* orbel_fget */
doublereal orbel_flon__(doublereal *e, doublereal *capn)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double sinh(doublereal);

    /* Local variables */
    static doublereal a, b, f;
    static integer i__;
    static doublereal x, a0, a1, b1, x2, fp, dx, sq, biga, bigb, diff;
    static integer iflag;

    /* Fortran I/O blocks */
    static cilist io___294 = { 0, 6, 0, 0, 0 };
    static cilist io___296 = { 0, 6, 0, 0, 0 };
    static cilist io___297 = { 0, 6, 0, 0, 0 };


/* ...  Inputs Only: */
/* ************************************************************************* */
/*                        SWIFT.INC */
/* ************************************************************************* */
/* Include file for SWIFT */

/* Author:  Hal Levison */
/* Date:    2/2/93 */
/* Last revision: 3/7/93 */
/* ...   Maximum array size */
/* you got it baby */
/* max number of planets, including th */
/* ...   Size of the test particle status flag */
/* max number of test particles */
/* ...   convergence criteria for danby */
/* ...    loop limits in the Laguerre attempts */
/* ...    A small number */
/* ...    trig stuff */
/* ------------------------------------------------------------------------- */
/* ...  Internals: */
/* ---- */
/* ...  Executable code */
/* Function to solve "Kepler's eqn" for F (here called */
/* x) for given e and CAPN. Only good for smallish CAPN */
    iflag = 0;
    if (*capn < 0.) {
	iflag = 1;
	*capn = -(*capn);
    }
    a1 = (1. - 1. / *e) * 6227020800.;
    a0 = *capn * -6227020800. / *e;
    b1 = a1;
/*  Set iflag nonzero if capn < 0., in which case solve for -capn */
/*  and change the sign of the final answer for F. */
/*  Begin with a reasonable guess based on solving the cubic for small F */
    a = (*e - 1.) * 6. / *e;
    b = *capn * -6. / *e;
    sq = sqrt(b * .25f * b + a * a * a / 27.);
    d__1 = b * -.5f + sq;
    biga = pow_dd(&d__1, &c_b299);
    d__1 = b * .5f + sq;
    bigb = -pow_dd(&d__1, &c_b299);
    x = biga + bigb;
/* 	write(6,*) 'cubic = ',x**3 +a*x +b */
    ret_val = x;
/* If capn is tiny (or zero) no need to go further than cubic even for */
/* e =1. */
    if (*capn < 4e-15) {
	goto L100;
    }
    for (i__ = 1; i__ <= 10; ++i__) {
	x2 = x * x;
	f = a0 + x * (a1 + x2 * (x2 * (x2 * (x2 * (x2 * (x2 + 156.) + 17160.) 
		+ 1235520.) + 51891840.) + 1037836800.));
	fp = b1 + x2 * (x2 * (x2 * (x2 * (x2 * (x2 * 13. + 1716.) + 154440.) 
		+ 8648640.) + 259459200.) + 3113510400.);
	dx = -f / fp;
/* 	  write(6,*) 'i,dx,x,f : ' */
/* 	  write(6,432) i,dx,x,f */
/* L432: */
	ret_val = x + dx;
/*   If we have converged here there's no point in going on */
	if (abs(dx) <= 4e-15) {
	    goto L100;
	}
	x = ret_val;
    }
/* Abnormal return here - we've gone thru the loop */
/* IMAX times without convergence */
    if (iflag == 1) {
	ret_val = -ret_val;
	*capn = -(*capn);
    }
    s_wsle(&io___294);
    do_lio(&c__9, &c__1, "FLON : RETURNING WITHOUT COMPLETE CONVERGENCE", (
	    ftnlen)45);
    e_wsle();
    diff = *e * sinh(ret_val) - ret_val - *capn;
    s_wsle(&io___296);
    do_lio(&c__9, &c__1, "N, F, ecc*sinh(F) - F - N : ", (ftnlen)28);
    e_wsle();
    s_wsle(&io___297);
    do_lio(&c__5, &c__1, (char *)&(*capn), (ftnlen)sizeof(doublereal));
    do_lio(&c__5, &c__1, (char *)&ret_val, (ftnlen)sizeof(doublereal));
    do_lio(&c__5, &c__1, (char *)&diff, (ftnlen)sizeof(doublereal));
    e_wsle();
    return ret_val;
/*  Normal return here, but check if capn was originally negative */
L100:
    if (iflag == 1) {
	ret_val = -ret_val;
	*capn = -(*capn);
    }
    return ret_val;
} /* orbel_flon__ */

/* ------------------------------------------------------------------ */

/* ********************************************************************** */
/*                    ORBEL_ZGET.F */
/* ********************************************************************** */
/*     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola */
/*          given Q (Fitz. notation.) */

/*             Input: */
/*                           q ==>  parabola mean anomaly. (real scalar) */
/*             Returns: */
/*                  orbel_zget ==>  eccentric anomaly. (real scalar) */

/*     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech." */
/*     REMARKS: For a parabola we can solve analytically. */
/*     AUTHOR: M. Duncan */
/*     DATE WRITTEN: May 11, 1992. */
/*     REVISIONS: May 27 - corrected it for negative Q and use power */
/* 	      series for small Q. */
/* ********************************************************************** */
/* orbel_flon */
doublereal orbel_zget__(doublereal *q)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal x, tmp;
    static integer iflag;

/* ...  Inputs Only: */
/* ************************************************************************* */
/*                        SWIFT.INC */
/* ************************************************************************* */
/* Include file for SWIFT */

/* Author:  Hal Levison */
/* Date:    2/2/93 */
/* Last revision: 3/7/93 */
/* ...   Maximum array size */
/* you got it baby */
/* max number of planets, including th */
/* ...   Size of the test particle status flag */
/* max number of test particles */
/* ...   convergence criteria for danby */
/* ...    loop limits in the Laguerre attempts */
/* ...    A small number */
/* ...    trig stuff */
/* ------------------------------------------------------------------------- */
/* ...  Internals: */
/* ---- */
/* ...  Executable code */
    iflag = 0;
    if (*q < 0.) {
	iflag = 1;
	*q = -(*q);
    }
    if (*q < .001) {
	ret_val = *q * (1. - *q * *q / 3. * (1. - *q * *q));
    } else {
/* Computing 2nd power */
	d__1 = *q;
	x = (*q * 3. + sqrt(d__1 * d__1 * 9. + 4.)) * .5;
	tmp = pow_dd(&x, &c_b317);
	ret_val = tmp - 1. / tmp;
    }
    if (iflag == 1) {
	ret_val = -ret_val;
	*q = -(*q);
    }
    return ret_val;
} /* orbel_zget__ */


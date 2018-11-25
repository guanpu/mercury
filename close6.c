/* close6.f -- translated by f2c (version 20181026).
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
static doublereal c_b83 = 0.;
static doublereal c_b84 = 11239424.;
static integer c__4 = 4;
static doublereal c_b114 = 3.141592653589793;
static doublereal c_b117 = 6.2831853071795862;
static doublereal c_b120 = 1.;
static doublereal c_b209 = 10.;
static integer c__8 = 8;
static integer c__7 = 7;
static doublereal c_b255 = 1461.;
static doublereal c_b256 = 365.25;
static doublereal c_b273 = .3333333333333333;
static integer c__5 = 5;
static doublereal c_b291 = .33333333333333331;

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      CLOSE6.FOR    (ErikSoft   5 June 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Makes output files containing details of close encounters that occurred */
/* during an integration using Mercury6 or higher. */

/* The user specifies the names of the required objects in the file close.in */

/* ------------------------------------------------------------------------------ */

/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;
    icilist ici__1;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer f_open(olist *), s_rsfe(cilist *), e_rsfe(void), f_clos(cllist *),
	     s_rsli(icilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsli(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_back(alist *);
    double d_lg10(doublereal *);
    integer s_wsfi(icilist *), e_wsfi(void), s_cmp(char *, char *, ftnlen, 
	    ftnlen);

    /* Local variables */
    extern /* Subroutine */ int mco_ov2x__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), m_formce__(integer *, char *, char *,
	     integer *, ftnlen, ftnlen), mio_jd_y__(doublereal *, integer *, 
	    integer *, doublereal *);
    static integer line_num__;
    static char c__[80*2000];
    static integer i__, j, k, l;
    static doublereal m[2000], a1, a2;
    static char c1[1];
    static doublereal e1, e2;
    static char master_id__[8*2000];
    static doublereal i1, i2, l1, n1;
    static integer firstflag;
    static doublereal p1, p2, n2, t0, t1, l2, v1[3], v2[3], x1[3];
    static integer precision;
    static doublereal x2[3], q1, q2;
    static integer timestyle;
    static char cc[80], id[8*2000];
    static doublereal gm, fr, fv;
    static char fin[5], mem[80*200];
    static doublereal phi;
    static integer lim[200]	/* was [2][100] */, master_unit__[2000];
    static doublereal rfac;
    static integer nbig;
    static doublereal jcen[3], dclo, mcen;
    static integer iclo, jclo, lmem[200];
    static doublereal rcen, time;
    static integer year;
    static doublereal vphi;
    static integer nsub;
    static doublereal rmax;
    static integer itmp, nsml;
    static char fout[250];
    static integer unit[2000];
    static logical test;
    static char type__[1], check[1];
    static integer nchar, algor, lenin;
    static doublereal theta;
    static integer nopen, nwait, month;
    static char style[1], header[250], infile__[250*50];
    static doublereal vtheta;
    static char string[250];
    static integer lenhead, allflag;
    extern /* Subroutine */ int mio_aei__(char *, char *, integer *, char *, 
	    integer *, char *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer ninfile;
    extern /* Subroutine */ int mio_err__(integer *, char *, integer *, char *
	    , integer *, char *, integer *, char *, integer *, ftnlen, ftnlen,
	     ftnlen, ftnlen), mio_spl__(integer *, char *, integer *, integer 
	    *, ftnlen);
    static integer nmaster;
    extern doublereal mio_c2fl__(char *, ftnlen), mio_c2re__(char *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    extern /* Subroutine */ int mco_x2el__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___4 = { 0, 14, 1, "(i3,1x,i2,1x,a80)", 0 };
    static cilist io___8 = { 0, 10, 0, "(a250)", 0 };
    static cilist io___13 = { 0, 10, 0, "(a250)", 0 };
    static cilist io___16 = { 0, 10, 0, "(a250)", 0 };
    static cilist io___24 = { 0, 10, 1, "(a250)", 0 };
    static cilist io___31 = { 1, 10, 1, "(3a1)", 0 };
    static cilist io___35 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___36 = { 0, 10, 0, "(3x,i2,a62,i1)", 0 };
    static cilist io___49 = { 1, 10, 0, "(a)", 0 };
    static icilist io___54 = { 0, fin+2, 0, "(i2)", 2, 1 };
    static cilist io___60 = { 1, 10, 0, "(3x,a70)", 0 };
    static cilist io___63 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___93 = { 0, 0, 0, fout, 0 };
    static cilist io___94 = { 0, 0, 0, fout, 0 };
    static cilist io___95 = { 0, 0, 0, fout, 0 };
    static cilist io___96 = { 0, 0, 0, fout, 0 };
    static cilist io___97 = { 0, 6, 0, "(2a,/,a,i10)", 0 };
    static cilist io___98 = { 0, 10, 1, "(a1)", 0 };



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
	s_wsfe(&io___3);
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
    i__1 = s_rsfe(&io___4);
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
    ioin__1.infilen = 8;
    ioin__1.infile = "close.in";
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
	o__1.ofnmlen = 8;
	o__1.ofnm = "close.in";
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    } else {
	mio_err__(&c__6, mem + 6400, &lmem[80], mem + 6960, &lmem[87], " ", &
		c__1, "close.in", &c__9, (ftnlen)80, (ftnlen)80, (ftnlen)1, (
		ftnlen)8);
    }

/* Read number of input files */
L30:
    s_rsfe(&io___8);
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
	s_rsfe(&io___13);
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

/* Read parameters used by this programme */
    timestyle = 1;
    for (j = 1; j <= 2; ++j) {
L50:
	s_rsfe(&io___16);
	do_fio(&c__1, string, (ftnlen)250);
	e_rsfe();
	if (*(unsigned char *)string == ')') {
	    goto L50;
	}
	mio_spl__(&c__250, string, &nsub, lim, (ftnlen)250);
	i__1 = lim[(nsub << 1) - 2] - 1;
	s_copy(c1, string + i__1, (ftnlen)1, lim[(nsub << 1) - 1] - i__1);
	if (j == 1 && (*(unsigned char *)c1 == 'd' || *(unsigned char *)c1 == 
		'D')) {
	    timestyle = 0;
	}
	if (j == 2 && (*(unsigned char *)c1 == 'y' || *(unsigned char *)c1 == 
		'Y')) {
	    timestyle += 2;
	}
    }

/* Read in the names of the objects for which orbital elements are required */
    nopen = 0;
    nwait = 0;
    nmaster = 0;
    m_formce__(&timestyle, fout, header, &lenhead, (ftnlen)250, (ftnlen)250);
L60:
    i__1 = s_rsfe(&io___24);
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
	mio_aei__(master_id__ + (nmaster - 1 << 3), ".clo", &master_unit__[
		nmaster - 1], header, &lenhead, mem, lmem, (ftnlen)8, (ftnlen)
		4, (ftnlen)250, (ftnlen)80);
    } else {
	++nwait;
	master_unit__[nmaster - 1] = -2;
    }
    goto L60;

L70:
/* If no objects are listed in CLOSE.IN assume that all objects are required */
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
	i__2 = s_rsfe(&io___31);
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
	    s_wsfe(&io___35);
	    do_fio(&c__1, " ERROR: This is an old style data file", (ftnlen)
		    38);
	    do_fio(&c__1, "        Try running m_close5.for instead.", (
		    ftnlen)41);
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
	    s_rsfe(&io___36);
	    do_fio(&c__1, (char *)&algor, (ftnlen)sizeof(integer));
	    do_fio(&c__1, cc, (ftnlen)62);
	    do_fio(&c__1, (char *)&precision, (ftnlen)sizeof(integer));
	    e_rsfe();

/* Decompress the time, number of objects, central mass and J components etc. */
	    time = mio_c2fl__(cc, (ftnlen)8);
	    if (firstflag == 0) {
		t0 = time;
		firstflag = 1;
	    }
	    nbig = (integer) (mio_c2re__(cc + 8, &c_b83, &c_b84, &c__3, (
		    ftnlen)8) + .5);
	    nsml = (integer) (mio_c2re__(cc + 11, &c_b83, &c_b84, &c__3, (
		    ftnlen)8) + .5);
	    mcen = mio_c2fl__(cc + 14, (ftnlen)8) * 2.959122082855911e-4;
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
		i__3 = s_rsfe(&io___49);
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
	    s_wsfi(&io___54);
	    do_fio(&c__1, (char *)&lenin, (ftnlen)sizeof(integer));
	    e_wsfi();

/* For each object decompress its name, code number, mass, spin and density */
	    i__2 = nbig + nsml;
	    for (j = 1; j <= i__2; ++j) {
		k = (integer) (mio_c2re__(c__ + (j - 1) * 80, &c_b83, &c_b84, 
			&c__3, (ftnlen)8) + .5);
		s_copy(id + (k - 1 << 3), c__ + ((j - 1) * 80 + 3), (ftnlen)8,
			 (ftnlen)8);
		m[k - 1] = mio_c2fl__(c__ + ((j - 1) * 80 + 11), (ftnlen)8) * 
			2.959122082855911e-4;

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
			    mio_aei__(master_id__ + (nmaster - 1 << 3), ".clo"
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

/*  IF  NORMAL  INPUT,  READ  COMPRESSED  DATA  ON  THE  CLOSE  ENCOUNTER */

	} else if (*(unsigned char *)type__ == 'b') {
	    ++line_num__;
	    i__2 = s_rsfe(&io___60);
	    if (i__2 != 0) {
		goto L666;
	    }
	    i__2 = do_fio(&c__1, cc, (ftnlen)70);
	    if (i__2 != 0) {
		goto L666;
	    }
	    i__2 = e_rsfe();
	    if (i__2 != 0) {
		goto L666;
	    }

/* Decompress time, distance and orbital variables for each object */
	    time = mio_c2fl__(cc, (ftnlen)8);
	    iclo = (integer) (mio_c2re__(cc + 8, &c_b83, &c_b84, &c__3, (
		    ftnlen)8) + .5);
	    jclo = (integer) (mio_c2re__(cc + 11, &c_b83, &c_b84, &c__3, (
		    ftnlen)8) + .5);
	    if (iclo > 2000 || jclo > 2000) {
		s_wsfe(&io___63);
		do_fio(&c__1, mem + 6400, lmem[80]);
		do_fio(&c__1, mem + 7120, lmem[89]);
		e_wsfe();
		s_stop("", (ftnlen)0);
	    }
	    dclo = mio_c2fl__(cc + 14, (ftnlen)8);
	    fr = mio_c2re__(cc + 22, &c_b83, &rfac, &c__4, (ftnlen)8);
	    theta = mio_c2re__(cc + 26, &c_b83, &c_b114, &c__4, (ftnlen)8);
	    phi = mio_c2re__(cc + 30, &c_b83, &c_b117, &c__4, (ftnlen)8);
	    fv = mio_c2re__(cc + 34, &c_b83, &c_b120, &c__4, (ftnlen)8);
	    vtheta = mio_c2re__(cc + 38, &c_b83, &c_b114, &c__4, (ftnlen)8);
	    vphi = mio_c2re__(cc + 42, &c_b83, &c_b117, &c__4, (ftnlen)8);
	    mco_ov2x__(&rcen, &rmax, &mcen, &m[iclo - 1], &fr, &theta, &phi, &
		    fv, &vtheta, &vphi, x1, &x1[1], &x1[2], v1, &v1[1], &v1[2]
		    );

	    fr = mio_c2re__(cc + 46, &c_b83, &rfac, &c__4, (ftnlen)8);
	    theta = mio_c2re__(cc + 50, &c_b83, &c_b114, &c__4, (ftnlen)8);
	    phi = mio_c2re__(cc + 54, &c_b83, &c_b117, &c__4, (ftnlen)8);
	    fv = mio_c2re__(cc + 58, &c_b83, &c_b120, &c__4, (ftnlen)8);
	    vtheta = mio_c2re__(cc + 62, &c_b83, &c_b114, &c__4, (ftnlen)8);
	    vphi = mio_c2re__(cc + 66, &c_b83, &c_b117, &c__4, (ftnlen)8);
	    mco_ov2x__(&rcen, &rmax, &mcen, &m[jclo - 1], &fr, &theta, &phi, &
		    fv, &vtheta, &vphi, x2, &x2[1], &x2[2], v2, &v2[1], &v2[2]
		    );

/* Convert to Keplerian elements */
	    gm = mcen + m[iclo - 1];
	    mco_x2el__(&gm, x1, &x1[1], &x1[2], v1, &v1[1], &v1[2], &q1, &e1, 
		    &i1, &p1, &n1, &l1);
	    a1 = q1 / (1. - e1);
	    gm = mcen + m[jclo - 1];
	    mco_x2el__(&gm, x2, &x2[1], &x2[2], v2, &v2[1], &v2[2], &q2, &e2, 
		    &i2, &p2, &n2, &l2);
	    a2 = q2 / (1. - e2);
	    i1 /= .017453292519943295;
	    i2 /= .017453292519943295;

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

/* Write encounter details to appropriate files */
	    if (timestyle == 1) {
		if (unit[iclo - 1] >= 10) {
		    io___93.ciunit = unit[iclo - 1];
		    s_wsfe(&io___93);
		    do_fio(&c__1, (char *)&year, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&month, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, id + (jclo - 1 << 3), (ftnlen)8);
		    do_fio(&c__1, (char *)&dclo, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&a1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&a2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}

		if (unit[jclo - 1] >= 10) {
		    io___94.ciunit = unit[jclo - 1];
		    s_wsfe(&io___94);
		    do_fio(&c__1, (char *)&year, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&month, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, id + (iclo - 1 << 3), (ftnlen)8);
		    do_fio(&c__1, (char *)&dclo, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&a2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&a1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
	    } else {
		if (unit[iclo - 1] >= 10) {
		    io___95.ciunit = unit[iclo - 1];
		    s_wsfe(&io___95);
		    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, id + (jclo - 1 << 3), (ftnlen)8);
		    do_fio(&c__1, (char *)&dclo, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&a1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&a2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(doublereal));
		    e_wsfe();
		}
		if (unit[jclo - 1] >= 10) {
		    io___96.ciunit = unit[jclo - 1];
		    s_wsfe(&io___96);
		    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, id + (iclo - 1 << 3), (ftnlen)8);
		    do_fio(&c__1, (char *)&dclo, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&a2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&a1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&e1, (ftnlen)sizeof(doublereal));
		    do_fio(&c__1, (char *)&i1, (ftnlen)sizeof(doublereal));
		    e_wsfe();
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
	s_wsfe(&io___97);
	do_fio(&c__1, mem + 9600, lmem[120]);
	do_fio(&c__1, infile__ + (i__ - 1) * 250, (ftnlen)60);
	do_fio(&c__1, mem + 8240, lmem[103]);
	do_fio(&c__1, (char *)&line_num__, (ftnlen)sizeof(integer));
	e_wsfe();
	*(unsigned char *)c1 = ' ';
	while(*(unsigned char *)c1 != 12) {
	    ++line_num__;
	    i__2 = s_rsfe(&io___98);
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

/* Move on to the next file containing close encounter data */
L900:
	cl__1.cerr = 0;
	cl__1.cunit = 10;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/* Close clo files */
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
		mio_aei__(master_id__ + (j - 1 << 3), ".clo", &master_unit__[
			j - 1], header, &lenhead, mem, lmem, (ftnlen)8, (
			ftnlen)4, (ftnlen)250, (ftnlen)80);
	    }
	}
	goto L90;
    }

    return 0;
} /* MAIN__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      M_FORMCE.FOR    (ErikSoft  30 November 1999) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */


/* ------------------------------------------------------------------------------ */

/* Subroutine */ int m_formce__(integer *timestyle, char *fout, char *header, 
	integer *lenhead, ftnlen fout_len, ftnlen header_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);



/* Input/Output */

/* ------------------------------------------------------------------------------ */

    if (*timestyle == 0 || *timestyle == 2) {
	s_copy(header, "    Time (days)    ", (ftnlen)19, (ftnlen)19);
	s_copy(header + 19, "  Object   dmin (AU)     a1       e1    ", (
		ftnlen)39, (ftnlen)40);
	s_copy(header + 58, "   i1       a2       e2       i2", (ftnlen)32, (
		ftnlen)32);
	*lenhead = 90;
	s_copy(fout, "(1x,f18.5,1x,a8,1x,f10.8,2(1x,f9.4,1x,f8.6,1x,f7.3))", (
		ftnlen)250, (ftnlen)52);
    } else {
	if (*timestyle == 1) {
	    s_copy(header, "     Year/Month/Day    ", (ftnlen)23, (ftnlen)23);
	    s_copy(header + 23, "  Object   dmin (AU)     a1       e1    ", (
		    ftnlen)39, (ftnlen)40);
	    s_copy(header + 62, "   i1       a2       e2       i2", (ftnlen)
		    32, (ftnlen)32);
	    *lenhead = 94;
	    s_copy(fout, "(1x,i10,1x,i2,1x,f8.5,1x,a8,1x,f10.8,", (ftnlen)37, 
		    (ftnlen)37);
	    s_copy(fout + 37, "2(1x,f9.4,1x,f8.6,1x,f7.3))", (ftnlen)27, (
		    ftnlen)27);
	} else {
	    s_copy(header, "    Time (years)   ", (ftnlen)19, (ftnlen)19);
	    s_copy(header + 19, "  Object   dmin (AU)     a1       e1    ", (
		    ftnlen)39, (ftnlen)40);
	    s_copy(header + 58, "   i1       a2       e2       i2", (ftnlen)
		    32, (ftnlen)32);
	    s_copy(fout, "(1x,f18.7,1x,a8,1x,f10.8,2(1x,f9.4,1x,f8.6,1x,f7.3"
		    "))", (ftnlen)250, (ftnlen)52);
	    *lenhead = 90;
	}
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* m_formce__ */


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

    r__ = *rcen * pow_dd(&c_b209, fr);
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
    static cilist io___161 = { 0, 6, 0, "(/,3a)", 0 };
    static cilist io___162 = { 0, 0, 0, "(/,30x,a8,//,a)", 0 };




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
	s_wsfe(&io___161);
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
	io___162.ciunit = *unitnum;
	s_wsfe(&io___162);
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

    x = mio_c2re__(c__, &c_b83, &c_b120, &c__7, (ftnlen)8);
    x = x * 2. - 1.;
    ex = (*(unsigned char *)&c__[7] + 256) % 256 - 144;
    d__1 = (doublereal) ex;
    ret_val = x * pow_dd(&c_b209, &d__1);

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
    static cilist io___167 = { 0, 6, 0, "(a)", 0 };
    static cilist io___168 = { 0, 0, 0, "(/,3a,/,2a)", 0 };




/* Input/Output */

/* ------------------------------------------------------------------------------ */

    s_wsfe(&io___167);
    do_fio(&c__1, " ERROR: Programme terminated.", (ftnlen)29);
    e_wsfe();
    io___168.ciunit = *unit;
    s_wsfe(&io___168);
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
		ce = d_sign(&c_b120, &ce);
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
	    cf = d_sign(&c_b120, &cf);
	}
	f = acos(cf);
	if (rv < 0.) {
	    f = 6.2831853071795862 - f;
	}
	*p = true__ - f;
	d__1 = *p + 6.2831853071795862 + 6.2831853071795862;
	*p = d_mod(&d__1, &c_b117);
    }

    if (*l < 0. && *e < 1.) {
	*l += 6.2831853071795862;
    }
    if (*l > 6.2831853071795862 && *e < 1.) {
	*l = d_mod(l, &c_b117);
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
    d__1 = d_mod(&x, &c_b255) + 1461.;
    y = d_int(&d__1);
    d__1 = d_mod(&y, &c_b256);
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
    static cilist io___236 = { 0, 6, 0, 0, 0 };


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
    s_wsle(&io___236);
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
    static cilist io___252 = { 0, 6, 0, 0, 0 };
    static cilist io___254 = { 0, 6, 0, 0, 0 };
    static cilist io___255 = { 0, 6, 0, 0, 0 };


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
    biga = pow_dd(&d__1, &c_b273);
    d__1 = b * .5f + sq;
    bigb = -pow_dd(&d__1, &c_b273);
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
    s_wsle(&io___252);
    do_lio(&c__9, &c__1, "FLON : RETURNING WITHOUT COMPLETE CONVERGENCE", (
	    ftnlen)45);
    e_wsle();
    diff = *e * sinh(ret_val) - ret_val - *capn;
    s_wsle(&io___254);
    do_lio(&c__9, &c__1, "N, F, ecc*sinh(F) - F - N : ", (ftnlen)28);
    e_wsle();
    s_wsle(&io___255);
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
	tmp = pow_dd(&x, &c_b291);
	ret_val = tmp - 1. / tmp;
    }
    if (iflag == 1) {
	ret_val = -ret_val;
	*q = -(*q);
    }
    return ret_val;
} /* orbel_zget__ */


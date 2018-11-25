/* mercury6_2.f -- translated by f2c (version 20181026).
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
static integer c__0 = 0;
static integer c__2 = 2;
static doublereal c_b58 = 6.2831853071795862;
static doublereal c_b95 = .3333333333333333;
static doublereal c_b98 = 0.;
static doublereal c_b99 = 11239423.99;
static doublereal c_b121 = 9.9e29;
static doublereal c_b150 = 1.;
static doublereal c_b167 = .1111111111111111;
static doublereal c_b176 = -2.15;
static doublereal c_b177 = -4.6142;
static doublereal c_b178 = 5.093;
static integer c__7 = 7;
static doublereal c_b186 = 10.;
static doublereal c_b196 = 3.141592653589793;
static integer c__9 = 9;
static integer c__5 = 5;
static integer c__3 = 3;
static integer c__6 = 6;
static integer c__8 = 8;
static integer c__150 = 150;
static integer c__80 = 80;
static integer c__23 = 23;
static doublereal c_b666 = .333333333333333;
static doublereal c_b957 = 1461.;
static doublereal c_b958 = 365.25;
static doublereal c_b1082 = .33333333333333331;

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MERCURY6_1.FOR    (ErikSoft   3 May 2002) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Mercury is a general-purpose N-body integration package for problems in */
/* celestial mechanics. */

/* ------------------------------------------------------------------------------ */
/* This package contains some subroutines taken from the Swift integration */
/* package by H.F.Levison and M.J.Duncan (1994) Icarus, vol 108, pp18. */
/* Routines taken from Swift have names beginning with `drift' or `orbel'. */

/* The standard symplectic (MVS) algorithm is described in J.Widsom and */
/* M.Holman (1991) Astronomical Journal, vol 102, pp1528. */

/* The hybrid symplectic algorithm is described in J.E.Chambers (1999) */
/* Monthly Notices of the RAS, vol 304, pp793. */

/* RADAU is described in E.Everhart (1985) in ``The Dynamics of Comets: */
/* Their Origin and Evolution'' p185-202, eds. A.Carusi & G.B.Valsecchi, */
/* pub. Reidel. */

/* The Bulirsch-Stoer algorithms are described in W.H.Press et al. (1992) */
/* ``Numerical Recipes in Fortran'', pub. Cambridge. */
/* ------------------------------------------------------------------------------ */

/* Variables: */
/* --------- */
/*  M      = mass (in solar masses) */
/*  XH     = coordinates (x,y,z) with respect to the central body (AU) */
/*  VH     = velocities (vx,vy,vz) with respect to the central body (AU/day) */
/*  S      = spin angular momentum (solar masses AU^2/day) */
/*  RHO    = physical density (g/cm^3) */
/*  RCEH   = close-encounter limit (Hill radii) */
/*  STAT   = status (0 => alive, <>0 => to be removed) */
/*  ID     = name of the object (8 characters) */
/*  CE     = close encounter status */
/*  NGF    = (1-3) cometary non-gravitational (jet) force parameters */
/*   "     =  (4)  beta parameter for radiation pressure and P-R drag */
/*  EPOCH  = epoch of orbit (days) */
/*  NBOD  = current number of bodies (INCLUDING the central object) */
/*  NBIG  =    "       "    " big bodies (ones that perturb everything else) */
/*  TIME  = current epoch (days) */
/*  TOUT  = time of next output evaluation */
/*  TDUMP = time of next data dump */
/*  TFUN  = time of next periodic effect (e.g. next check for ejections) */
/*  H     = current integration timestep (days) */
/*  EN(1) = initial energy of the system */
/*  " (2) = current    "    "  "    " */
/*  " (3) = energy change due to collisions, ejections etc. */
/*  AM(1,2,3) = as above but for angular momentum */

/* Integration Parameters : */
/* ---------------------- */
/*  ALGOR = 1  ->  Mixed-variable symplectic */
/*          2  ->  Bulirsch-Stoer integrator */
/*          3  ->         "           "      (conservative systems only) */
/*          4  ->  RA15 `radau' integrator */
/*          10 ->  Hybrid MVS/BS (democratic-heliocentric coords) */
/*          11 ->  Close-binary hybrid (close-binary coords) */
/*          12 ->  Wide-binary hybrid (wide-binary coords) */

/* TSTART = epoch of first required output (days) */
/* TSTOP  =   "      final required output ( "  ) */
/* DTOUT  = data output interval           ( "  ) */
/* DTDUMP = data-dump interval             ( "  ) */
/* DTFUN  = interval for other periodic effects (e.g. check for ejections) */
/*  H0    = initial integration timestep (days) */
/*  TOL   = Integrator tolerance parameter (approx. error per timestep) */
/*  RMAX  = heliocentric distance at which objects are considered ejected (AU) */
/*  RCEN  = radius of central body (AU) */
/*  JCEN(1,2,3) = J2,J4,J6 for central body (units of RCEN^i for Ji) */

/* Options: */
/*  OPT(1) = close-encounter option (0=stop after an encounter, 1=continue) */
/*  OPT(2) = collision option (0=no collisions, 1=merge, 2=merge+fragment) */
/*  OPT(3) = time style (0=days 1=Greg.date 2/3=days/years w/respect to start) */
/*  OPT(4) = o/p precision (1,2,3 = 4,9,15 significant figures) */
/*  OPT(5) = < Not used at present > */
/*  OPT(6) = < Not used at present > */
/*  OPT(7) = apply post-Newtonian correction? (0=no, 1=yes) */
/*  OPT(8) = apply user-defined force routine mfo_user? (0=no, 1=yes) */

/* File variables : */
/* -------------- */
/*  OUTFILE  (1) = osculating coordinates/velocities and masses */
/*     "     (2) = close encounter details */
/*     "     (3) = information file */
/*  DUMPFILE (1) = Big-body data */
/*     "     (2) = Small-body data */
/*     "     (3) = integration parameters */
/*     "     (4) = restart file */

/* Flags : */
/* ----- */
/*  NGFLAG = do any bodies experience non-grav. forces? */
/*                            ( 0 = no non-grav forces) */
/*                              1 = cometary jets only */
/*                              2 = radiation pressure/P-R drag only */
/*                              3 = both */
/*  OPFLAG = integration mode (-2 = synchronising epochs) */
/*                             -1 = integrating towards start epoch */
/*                              0 = main integration, normal output */
/*                              1 = main integration, full output */

/* ------------------------------------------------------------------------------ */

/* Main program */ int MAIN__(void)
{
    /* Initialized data */

    static integer opt[8] = { 0,1,1,2,0,1,0,0 };

    /* Format strings */
    static char fmt_231[] = "(/,a,1p1e12.5)";
    static char fmt_232[] = "(a,1p1e12.5)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void), f_clos(cllist *);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int mco_iden__();
    extern /* Subroutine */ int mal_hcon__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, char *, doublereal *, 
	    integer *, integer *, integer *, char *, char *, char *, integer *
	    , U_fp, U_fp, U_fp, ftnlen, ftnlen, ftnlen, ftnlen), mal_hvar__(
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     char *, doublereal *, integer *, integer *, integer *, char *, 
	    char *, char *, integer *, U_fp, ftnlen, ftnlen, ftnlen, ftnlen);
    static char dumpfile[80*4];
    extern /* Subroutine */ int mio_dump__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, char *, doublereal *, 
	    doublereal *, integer *, integer *, char *, char *, integer *, 
	    ftnlen, ftnlen, ftnlen), mxx_sync__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, char *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen);
    static integer j;
    extern /* Subroutine */ int mco_h2mvs__();
    static doublereal m[2000];
    extern /* Subroutine */ int mco_mvs2h__();
    static doublereal s[6000]	/* was [3][2000] */, h0;
    static char id[8*2000];
    static doublereal am[3], en[3], vh[6000]	/* was [3][2000] */, xh[6000]	
	    /* was [3][2000] */, ngf[8000]	/* was [4][2000] */;
    static char mem[80*200];
    static doublereal rho[2000], tol;
    static integer nbig;
    static doublereal jcen[3], rceh[2000];
    static integer nbod;
    static doublereal rcen;
    static integer lmem[200];
    static doublereal time;
    static integer nfun;
    static doublereal rmax;
    static integer stat[2000];
    static doublereal cefac, epoch[2000];
    static integer algor, ndump;
    static doublereal dtout, tstop;
    static integer ngflag, opflag;
    extern /* Subroutine */ int mio_in__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, char *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, char *, char *, 
	    integer *, char *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int mdt_hy__();
    extern /* Subroutine */ int mxx_en__(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal tstart;
    extern /* Subroutine */ int mdt_bs1__(), mdt_bs2__(), mdt_ra15__();
    static char outfile[80*3];
    extern /* Subroutine */ int mdt_mvs__(), mco_dh2h__(), mco_h2dh__();

    /* Fortran I/O blocks */
    static cilist io___35 = { 0, 23, 0, "(/,a)", 0 };
    static cilist io___36 = { 0, 6, 0, "(a)", 0 };
    static cilist io___37 = { 0, 23, 0, "(/,a,/)", 0 };
    static cilist io___38 = { 0, 6, 0, "(a)", 0 };
    static cilist io___40 = { 0, 23, 0, "(/,a)", 0 };
    static cilist io___41 = { 0, 23, 0, fmt_231, 0 };
    static cilist io___42 = { 0, 23, 0, fmt_232, 0 };
    static cilist io___43 = { 0, 23, 0, fmt_231, 0 };
    static cilist io___44 = { 0, 23, 0, fmt_232, 0 };
    static cilist io___45 = { 0, 6, 0, "(a)", 0 };



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

/* Get initial conditions and integration parameters */
    mio_in__(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, &rmax, &rcen, 
	    jcen, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s, 
	    rho, rceh, stat, id, epoch, ngf, opt, &opflag, &ngflag, outfile, 
	    dumpfile, lmem, mem, (ftnlen)8, (ftnlen)80, (ftnlen)80, (ftnlen)
	    80);

/* If this is a new integration, integrate all the objects to a common epoch. */
    if (opflag == -2) {
L20:
	o__1.oerr = 1;
	o__1.ounit = 23;
	o__1.ofnmlen = 80;
	o__1.ofnm = outfile + 160;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = "append";
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L20;
	}
	s_wsfe(&io___35);
	do_fio(&c__1, mem + 4320, lmem[54]);
	e_wsfe();
	s_wsfe(&io___36);
	do_fio(&c__1, mem + 4320, lmem[54]);
	e_wsfe();
	mxx_sync__(&time, &tstart, &h0, &tol, jcen, &nbod, &nbig, m, xh, vh, 
		s, rho, rceh, stat, id, epoch, ngf, opt, &ngflag, (ftnlen)8);
	s_wsfe(&io___37);
	do_fio(&c__1, mem + 4400, lmem[55]);
	e_wsfe();
	s_wsfe(&io___38);
	do_fio(&c__1, mem + 4400, lmem[55]);
	e_wsfe();
	opflag = -1;
	cl__1.cerr = 0;
	cl__1.cunit = 23;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/* Main integration */
    if (algor == 1) {
	mal_hcon__(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &
		rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, 
		xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, 
		outfile, dumpfile, mem, lmem, (U_fp)mdt_mvs__, (U_fp)
		mco_h2mvs__, (U_fp)mco_mvs2h__, (ftnlen)8, (ftnlen)80, (
		ftnlen)80, (ftnlen)80);
    }

    if (algor == 9) {
	mal_hcon__(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &
		rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, 
		xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, 
		outfile, dumpfile, mem, lmem, (U_fp)mdt_mvs__, (U_fp)
		mco_iden__, (U_fp)mco_iden__, (ftnlen)8, (ftnlen)80, (ftnlen)
		80, (ftnlen)80);
    }

    if (algor == 2) {
	mal_hvar__(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &
		rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, 
		xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, 
		outfile, dumpfile, mem, lmem, (U_fp)mdt_bs1__, (ftnlen)8, (
		ftnlen)80, (ftnlen)80, (ftnlen)80);
    }

    if (algor == 3) {
	mal_hvar__(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &
		rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, 
		xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, 
		outfile, dumpfile, mem, lmem, (U_fp)mdt_bs2__, (ftnlen)8, (
		ftnlen)80, (ftnlen)80, (ftnlen)80);
    }

    if (algor == 4) {
	mal_hvar__(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &
		rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, 
		xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, 
		outfile, dumpfile, mem, lmem, (U_fp)mdt_ra15__, (ftnlen)8, (
		ftnlen)80, (ftnlen)80, (ftnlen)80);
    }

    if (algor == 10) {
	mal_hcon__(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &
		rcen, &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, 
		xh, vh, s, rho, rceh, stat, id, ngf, opt, &opflag, &ngflag, 
		outfile, dumpfile, mem, lmem, (U_fp)mdt_hy__, (U_fp)
		mco_h2dh__, (U_fp)mco_dh2h__, (ftnlen)8, (ftnlen)80, (ftnlen)
		80, (ftnlen)80);
    }

/* Do a final data dump */
    i__1 = nbod;
    for (j = 2; j <= i__1; ++j) {
	epoch[j - 1] = time;
    }
    mio_dump__(&time, &tstart, &tstop, &dtout, &algor, &h0, &tol, jcen, &rcen,
	     &rmax, en, am, &cefac, &ndump, &nfun, &nbod, &nbig, m, xh, vh, s,
	     rho, rceh, stat, id, ngf, epoch, opt, &opflag, dumpfile, mem, 
	    lmem, (ftnlen)8, (ftnlen)80, (ftnlen)80);

/* Calculate and record the overall change in energy and ang. momentum */
L50:
    o__1.oerr = 1;
    o__1.ounit = 23;
    o__1.ofnmlen = 80;
    o__1.ofnm = outfile + 160;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = "append";
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L50;
    }
    s_wsfe(&io___40);
    do_fio(&c__1, mem + 4480, lmem[56]);
    e_wsfe();
    mxx_en__(jcen, &nbod, &nbig, m, xh, vh, s, &en[1], &am[1]);

    s_wsfe(&io___41);
    do_fio(&c__1, mem + 4560, lmem[57]);
    d__2 = (d__1 = (en[1] + en[2] - en[0]) / en[0], abs(d__1));
    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___42);
    do_fio(&c__1, mem + 4640, lmem[58]);
    d__2 = (d__1 = (am[1] + am[2] - am[0]) / am[0], abs(d__1));
    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___43);
    do_fio(&c__1, mem + 4720, lmem[59]);
    d__2 = (d__1 = en[2] / en[0], abs(d__1));
    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___44);
    do_fio(&c__1, mem + 4800, lmem[60]);
    d__2 = (d__1 = am[2] / am[0], abs(d__1));
    do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 23;
    cl__1.csta = 0;
    f_clos(&cl__1);
    s_wsfe(&io___45);
    do_fio(&c__1, mem + 4480, lmem[56]);
    e_wsfe();

/* ------------------------------------------------------------------------------ */

    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MFO_USER.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Applies an arbitrary force, defined by the user. */

/* If using with the symplectic algorithm MAL_MVS, the force should be */
/* small compared with the force from the central object. */
/* If using with the conservative Bulirsch-Stoer algorithm MAL_BS2, the */
/* force should not be a function of the velocities. */

/* N.B. All coordinates and velocities must be with respect to central body */
/* === */
/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mfo_user__(doublereal *time, doublereal *jcen, integer *
	nbod, integer *nbig, doublereal *m, doublereal *x, doublereal *v, 
	doublereal *a)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;



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
    --jcen;
    a -= 4;
    v -= 4;
    x -= 4;
    --m;

    /* Function Body */
    i__1 = *nbod;
    for (j = 1; j <= i__1; ++j) {
	a[j * 3 + 1] = 0.;
	a[j * 3 + 2] = 0.;
	a[j * 3 + 3] = 0.;
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mfo_user__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MAL_HVAR.FOR    (ErikSoft   4 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Does an integration using a variable-timestep integration algorithm. The */
/* particular integrator routine is ONESTEP and the algorithm must use */
/* coordinates with respect to the central body. */

/* N.B. This routine is also called by the synchronisation routine mxx_sync, */
/* ===  in which case OPFLAG = -2. Beware when making changes involving OPFLAG. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mal_hvar__(doublereal *time, doublereal *tstart, 
	doublereal *tstop, doublereal *dtout, integer *algor, doublereal *h0, 
	doublereal *tol, doublereal *jcen, doublereal *rcen, doublereal *rmax,
	 doublereal *en, doublereal *am, doublereal *cefac, integer *ndump, 
	integer *nfun, integer *nbod, integer *nbig, doublereal *m, 
	doublereal *xh, doublereal *vh, doublereal *s, doublereal *rho, 
	doublereal *rceh, integer *stat, char *id, doublereal *ngf, integer *
	opt, integer *opflag, integer *ngflag, char *outfile, char *dumpfile, 
	char *mem, integer *lmem, S_fp onestep, ftnlen id_len, ftnlen 
	outfile_len, ftnlen dumpfile_len, ftnlen mem_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    extern /* Subroutine */ int mce_hill__(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), 
	    mco_iden__(doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    mce_coll__(doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     char *, integer *, char *, integer *, char *, ftnlen, ftnlen, 
	    ftnlen), mce_cent__(doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), mce_init__(doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, char *, integer *, ftnlen, 
	    ftnlen), mce_stat__(doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    char *, char *, integer *, ftnlen, ftnlen), mxx_ejec__(doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, char *, 
	    integer *, integer *, char *, char *, integer *, ftnlen, ftnlen, 
	    ftnlen), mio_dump__(doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, char *, doublereal *, 
	    doublereal *, integer *, integer *, char *, char *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer stopflag;
    extern /* Subroutine */ int mxx_elim__(integer *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, char *, char 
	    *, integer *, char *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal a[2000], h__;
    static integer i__, j, k, n;
    static doublereal v0[6000]	/* was [3][2000] */, x0[6000]	/* was [3][
	    2000] */;
    static integer ce[2000], ice[2000], jce[2000], nce;
    static doublereal rce[2000], tmp0, hdid, dclo[50];
    static integer iclo[50], chit[50], nclo, ihit[50], jhit[50], jclo[50];
    static doublereal dhit[50];
    static integer nhit;
    static doublereal tclo[50], tlog, thit[50];
    static integer itmp;
    static doublereal tfun, tout, thit1, epoch[2000], dtfun, rcrit[2000], 
	    tdump, rphys[2000];
    static integer ejflag;
    extern /* Subroutine */ int mio_ce__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, char *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     char *, integer *, char *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer dtflag;
    static doublereal tsmall, dtdump;
    extern /* Subroutine */ int mxx_en__(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal ixvclo[300]	/* was [6][50] */, jxvclo[300]	/* 
	    was [6][50] */;
    extern /* Subroutine */ int mfo_all__();
    extern /* Subroutine */ int mio_log__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, char *, integer *, ftnlen);
    static integer nowflag;
    extern /* Subroutine */ int mio_out__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     char *, integer *, integer *, integer *, char *, ftnlen, ftnlen);
    static integer nstored;



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

/* Initialize variables. DTFLAG = 0 implies first ever call to ONESTEP */
    /* Parameter adjustments */
    --jcen;
    --en;
    --am;
    ngf -= 5;
    id -= 8;
    --stat;
    --rceh;
    --rho;
    s -= 4;
    vh -= 4;
    xh -= 4;
    --m;
    --opt;
    outfile -= 80;
    dumpfile -= 80;
    mem -= 80;
    --lmem;

    /* Function Body */
    *dtout = abs(*dtout);
    dtdump = abs(*h0) * *ndump;
    dtfun = abs(*h0) * *nfun;
    dtflag = 0;
    nstored = 0;
    tsmall = *h0 * 1e-8;
    h__ = *h0;
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	ce[j - 1] = 0;
    }

/* Calculate close-encounter limits and physical radii for massive bodies */
    mce_init__(tstart, algor, h0, &jcen[1], rcen, rmax, cefac, nbod, nbig, &m[
	    1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1], rphys, rce, rcrit, 
	    id + 8, &opt[1], outfile + 160, &c__1, (ftnlen)8, (ftnlen)80);

/* Set up time of next output, times of previous dump, log and periodic effect */
    if (*opflag == -1) {
	tout = *tstart;
    } else {
	n = (integer) ((d__1 = *time - *tstart, abs(d__1)) / *dtout) + 1;
	d__1 = (doublereal) n;
	d__2 = *tstop - *tstart;
	tout = *tstart + *dtout * d_sign(&d__1, &d__2);
	if ((*tstop - *tstart) * (tout - *tstop) > 0.) {
	    tout = *tstop;
	}
    }
    tdump = *time;
    tfun = *time;
    tlog = *time;

/* ------------------------------------------------------------------------------ */

/*  MAIN  LOOP  STARTS  HERE */

L100:

/* Is it time for output ? */
    if ((d__1 = tout - *time, abs(d__1)) < abs(tsmall) && *opflag >= -1) {

/* Beware: the integration may change direction at this point!!!! */
	if (*opflag == -1) {
	    dtflag = 0;
	}

/* Output data for all bodies */
	mio_out__(time, &jcen[1], rcen, rmax, nbod, nbig, &m[1], &xh[4], &vh[
		4], &s[4], &rho[1], &stat[1], id + 8, &opt[1], opflag, algor, 
		outfile + 80, (ftnlen)8, (ftnlen)80);
	mio_ce__(time, tstart, rcen, rmax, nbod, nbig, &m[1], &stat[1], id + 
		8, &c__0, iclo, jclo, &opt[1], &stopflag, tclo, dclo, ixvclo, 
		jxvclo, mem + 80, &lmem[1], outfile + 80, &nstored, &c__0, (
		ftnlen)8, (ftnlen)80, (ftnlen)80);
	tmp0 = *tstop - tout;
/* Computing MIN */
	d__2 = abs(tmp0), d__3 = abs(*dtout);
	d__1 = min(d__2,d__3);
	tout += d_sign(&d__1, &tmp0);

/* Update the data dump files */
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    epoch[j - 1] = *time;
	}
	mio_dump__(time, tstart, tstop, dtout, algor, &h__, tol, &jcen[1], 
		rcen, rmax, &en[1], &am[1], cefac, ndump, nfun, nbod, nbig, &
		m[1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1], &stat[1], id 
		+ 8, &ngf[5], epoch, &opt[1], opflag, dumpfile + 80, mem + 80,
		 &lmem[1], (ftnlen)8, (ftnlen)80, (ftnlen)80);
	tdump = *time;
    }

/* If integration has finished return to the main part of programme */
    if ((d__1 = *tstop - *time, abs(d__1)) <= abs(tsmall) && *opflag != -1) {
	return 0;
    }

/* Set the timestep */
    if (*opflag == -1) {
	tmp0 = *tstart - *time;
    }
    if (*opflag == -2) {
	tmp0 = *tstop - *time;
    }
    if (*opflag >= 0) {
	tmp0 = tout - *time;
    }
/* Computing MAX */
/* Computing MIN */
    d__3 = abs(tmp0), d__4 = abs(h__);
    d__2 = min(d__3,d__4);
    d__1 = max(d__2,tsmall);
    h__ = d_sign(&d__1, &tmp0);

/* Save the current coordinates and velocities */
    mco_iden__(time, &jcen[1], nbod, nbig, &h__, &m[1], &xh[4], &vh[4], x0, 
	    v0, &ngf[5], ngflag, &opt[1]);

/* Advance one timestep */
    (*onestep)(time, &h__, &hdid, tol, &jcen[1], nbod, nbig, &m[1], &xh[4], &
	    vh[4], &s[4], rphys, rcrit, &ngf[5], &stat[1], &dtflag, ngflag, &
	    opt[1], &nce, ice, jce, (U_fp)mfo_all__);
    *time += hdid;

/* Check if close encounters or collisions occurred */
    nclo = 0;
    mce_stat__(time, &h__, rcen, nbod, nbig, &m[1], x0, v0, &xh[4], &vh[4], 
	    rce, rphys, &nclo, iclo, jclo, dclo, tclo, ixvclo, jxvclo, &nhit, 
	    ihit, jhit, chit, dhit, thit, &thit1, &nowflag, &stat[1], outfile 
	    + 240, mem + 80, &lmem[1], (ftnlen)80, (ftnlen)80);

/* ------------------------------------------------------------------------------ */

/*  CLOSE  ENCOUNTERS */

/* If encounter minima occurred, output details and decide whether to stop */
    if (nclo > 0 && *opflag >= -1) {
	itmp = 1;
	if (nhit != 0) {
	    itmp = 0;
	}
	mio_ce__(time, tstart, rcen, rmax, nbod, nbig, &m[1], &stat[1], id + 
		8, &nclo, iclo, jclo, &opt[1], &stopflag, tclo, dclo, ixvclo, 
		jxvclo, mem + 80, &lmem[1], outfile + 80, &nstored, &itmp, (
		ftnlen)8, (ftnlen)80, (ftnlen)80);
	if (stopflag == 1) {
	    return 0;
	}
    }

/* ------------------------------------------------------------------------------ */

/*  COLLISIONS */

/* If a collision occurred, output details and resolve the collision */
    if (nhit > 0 && opt[2] != 0) {
	i__1 = nhit;
	for (k = 1; k <= i__1; ++k) {
	    if (chit[k - 1] == 1) {
		i__ = ihit[k - 1];
		j = jhit[k - 1];
		mce_coll__(&thit[k - 1], tstart, &en[3], &jcen[1], &i__, &j, 
			nbod, nbig, &m[1], &xh[4], &vh[4], &s[4], rphys, &
			stat[1], id + 8, &opt[1], mem + 80, &lmem[1], outfile 
			+ 240, (ftnlen)8, (ftnlen)80, (ftnlen)80);
	    }
	}

/* Remove lost objects, reset flags and recompute Hill and physical radii */
	mxx_elim__(nbod, nbig, &m[1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1]
		, rcrit, &ngf[5], &stat[1], id + 8, mem + 80, &lmem[1], 
		outfile + 240, &itmp, (ftnlen)8, (ftnlen)80, (ftnlen)80);
	dtflag = 1;
	if (*opflag >= 0) {
	    *opflag = 1;
	}
	mce_init__(tstart, algor, h0, &jcen[1], rcen, rmax, cefac, nbod, nbig,
		 &m[1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1], rphys, rce, 
		rcrit, id + 8, &opt[1], outfile + 160, &c__1, (ftnlen)8, (
		ftnlen)80);
    }

/* ------------------------------------------------------------------------------ */

/*  COLLISIONS  WITH  CENTRAL  BODY */

/* Check for collisions */
    mce_cent__(time, &hdid, rcen, &jcen[1], &c__2, nbod, nbig, &m[1], x0, v0, 
	    &xh[4], &vh[4], &nhit, jhit, thit, dhit, algor, &ngf[5], ngflag);

/* Resolve the collisions */
    if (nhit > 0) {
	i__1 = nhit;
	for (k = 1; k <= i__1; ++k) {
	    i__ = 1;
	    j = jhit[k - 1];
	    mce_coll__(&thit[k - 1], tstart, &en[3], &jcen[1], &i__, &j, nbod,
		     nbig, &m[1], &xh[4], &vh[4], &s[4], rphys, &stat[1], id 
		    + 8, &opt[1], mem + 80, &lmem[1], outfile + 240, (ftnlen)
		    8, (ftnlen)80, (ftnlen)80);
	}

/* Remove lost objects, reset flags and recompute Hill and physical radii */
	mxx_elim__(nbod, nbig, &m[1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1]
		, rcrit, &ngf[5], &stat[1], id + 8, mem + 80, &lmem[1], 
		outfile + 240, &itmp, (ftnlen)8, (ftnlen)80, (ftnlen)80);
	dtflag = 1;
	if (*opflag >= 0) {
	    *opflag = 1;
	}
	mce_init__(tstart, algor, h0, &jcen[1], rcen, rmax, cefac, nbod, nbig,
		 &m[1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1], rphys, rce, 
		rcrit, id + 8, &opt[1], outfile + 160, &c__0, (ftnlen)8, (
		ftnlen)80);
    }

/* ------------------------------------------------------------------------------ */

/*  DATA  DUMP  AND  PROGRESS  REPORT */

/* Do the data dump */
    if ((d__1 = *time - tdump, abs(d__1)) >= abs(dtdump) && *opflag >= -1) {
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    epoch[j - 1] = *time;
	}
	mio_ce__(time, tstart, rcen, rmax, nbod, nbig, &m[1], &stat[1], id + 
		8, &c__0, iclo, jclo, &opt[1], &stopflag, tclo, dclo, ixvclo, 
		jxvclo, mem + 80, &lmem[1], outfile + 80, &nstored, &c__0, (
		ftnlen)8, (ftnlen)80, (ftnlen)80);
	mio_dump__(time, tstart, tstop, dtout, algor, &h__, tol, &jcen[1], 
		rcen, rmax, &en[1], &am[1], cefac, ndump, nfun, nbod, nbig, &
		m[1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1], &stat[1], id 
		+ 8, &ngf[5], epoch, &opt[1], opflag, dumpfile + 80, mem + 80,
		 &lmem[1], (ftnlen)8, (ftnlen)80, (ftnlen)80);
	tdump = *time;
    }

/* Write a progress report to the log file */
    if ((d__1 = *time - tlog, abs(d__1)) >= abs(dtdump) && *opflag >= 0) {
	mxx_en__(&jcen[1], nbod, nbig, &m[1], &xh[4], &vh[4], &s[4], &en[2], &
		am[2]);
	mio_log__(time, tstart, &en[1], &am[1], &opt[1], mem + 80, &lmem[1], (
		ftnlen)80);
	tlog = *time;
    }

/* ------------------------------------------------------------------------------ */

/*  CHECK  FOR  EJECTIONS  AND  DO  OTHER  PERIODIC  EFFECTS */

    if ((d__1 = *time - tfun, abs(d__1)) >= abs(dtfun) && *opflag >= -1) {

/* Recompute close encounter limits, to allow for changes in Hill radii */
	mce_hill__(nbod, &m[1], &xh[4], &vh[4], rce, a);
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    rce[j - 1] *= rceh[j];
	}

/* Check for ejections */
	mxx_ejec__(time, tstart, rmax, &en[1], &am[1], &jcen[1], &c__2, nbod, 
		nbig, &m[1], &xh[4], &vh[4], &s[4], &stat[1], id + 8, &opt[1],
		 &ejflag, outfile + 240, mem + 80, &lmem[1], (ftnlen)8, (
		ftnlen)80, (ftnlen)80);

/* Remove lost objects, reset flags and recompute Hill and physical radii */
	if (ejflag != 0) {
	    mxx_elim__(nbod, nbig, &m[1], &xh[4], &vh[4], &s[4], &rho[1], &
		    rceh[1], rcrit, &ngf[5], &stat[1], id + 8, mem + 80, &
		    lmem[1], outfile + 240, &itmp, (ftnlen)8, (ftnlen)80, (
		    ftnlen)80);
	    dtflag = 1;
	    if (*opflag >= 0) {
		*opflag = 1;
	    }
	    mce_init__(tstart, algor, h0, &jcen[1], rcen, rmax, cefac, nbod, 
		    nbig, &m[1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1], 
		    rphys, rce, rcrit, id + 8, &opt[1], outfile + 160, &c__0, 
		    (ftnlen)8, (ftnlen)80);
	}
	tfun = *time;
    }

/* Go on to the next time step */
    goto L100;

/* ------------------------------------------------------------------------------ */

} /* mal_hvar__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MAL_HCON.FOR    (ErikSoft   28 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Does an integration using an integrator with a constant stepsize H. */
/* Input and output to this routine use coordinates XH, and velocities VH, */
/* with respect to the central body, but the integration algorithm uses */
/* its own internal coordinates X, and velocities V. */

/* The programme uses the transformation routines COORD and BCOORD to change */
/* to and from the internal coordinates, respectively. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mal_hcon__(doublereal *time, doublereal *tstart, 
	doublereal *tstop, doublereal *dtout, integer *algor, doublereal *h0, 
	doublereal *tol, doublereal *jcen, doublereal *rcen, doublereal *rmax,
	 doublereal *en, doublereal *am, doublereal *cefac, integer *ndump, 
	integer *nfun, integer *nbod, integer *nbig, doublereal *m, 
	doublereal *xh, doublereal *vh, doublereal *s, doublereal *rho, 
	doublereal *rceh, integer *stat, char *id, doublereal *ngf, integer *
	opt, integer *opflag, integer *ngflag, char *outfile, char *dumpfile, 
	char *mem, integer *lmem, S_fp onestep, S_fp coord, S_fp bcoord, 
	ftnlen id_len, ftnlen outfile_len, ftnlen dumpfile_len, ftnlen 
	mem_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    extern /* Subroutine */ int mce_hill__(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), 
	    mco_iden__(doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    mce_cent__(doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *), 
	    mce_coll__(doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     char *, integer *, char *, integer *, char *, ftnlen, ftnlen, 
	    ftnlen), mce_init__(doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, char *, integer *, char *, integer *, ftnlen, 
	    ftnlen), mxx_ejec__(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, char *, integer *, integer *, char *, char *, integer 
	    *, ftnlen, ftnlen, ftnlen), mio_dump__(doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, char *, doublereal *, 
	    doublereal *, integer *, integer *, char *, char *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer stopflag;
    extern /* Subroutine */ int mxx_elim__(integer *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, char *, char 
	    *, integer *, char *, integer *, ftnlen, ftnlen, ftnlen);
    static doublereal a[2000];
    static integer i__, j, k, n;
    static doublereal v[6000]	/* was [3][2000] */, x[6000]	/* was [3][
	    2000] */, vh0[6000]	/* was [3][2000] */, xh0[6000]	/* was [3][
	    2000] */, rce[2000], hby2, tmp0, dclo[50];
    static integer iclo[50], nclo, jclo[50];
    static doublereal dhit[50];
    static integer jhit[50];
    static doublereal tclo[50];
    static integer nhit;
    static doublereal tlog, thit[50];
    static integer itmp;
    static doublereal tfun, tout, epoch[2000], dtfun, rcrit[2000], tdump, 
	    rphys[2000];
    static integer ejflag;
    extern /* Subroutine */ int mio_ce__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, char *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     char *, integer *, char *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer dtflag;
    static doublereal dtdump;
    extern /* Subroutine */ int mxx_en__(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal ixvclo[300]	/* was [6][50] */, jxvclo[300]	/* 
	    was [6][50] */;
    static integer colflag;
    extern /* Subroutine */ int mio_log__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, char *, integer *, ftnlen),
	     mio_out__(doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, char *, integer *, integer 
	    *, integer *, char *, ftnlen, ftnlen);
    static integer nstored;



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

/* Initialize variables. DTFLAG = 0/2: first call ever/normal */
    /* Parameter adjustments */
    --jcen;
    --en;
    --am;
    ngf -= 5;
    id -= 8;
    --stat;
    --rceh;
    --rho;
    s -= 4;
    vh -= 4;
    xh -= 4;
    --m;
    --opt;
    outfile -= 80;
    dumpfile -= 80;
    mem -= 80;
    --lmem;

    /* Function Body */
    *dtout = abs(*dtout);
    dtdump = abs(*h0) * *ndump;
    dtfun = abs(*h0) * *nfun;
    dtflag = 0;
    nstored = 0;
    hby2 = abs(*h0) * .500001;

/* Calculate close-encounter limits and physical radii */
    mce_init__(tstart, algor, h0, &jcen[1], rcen, rmax, cefac, nbod, nbig, &m[
	    1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1], rphys, rce, rcrit, 
	    id + 8, &opt[1], outfile + 160, &c__1, (ftnlen)8, (ftnlen)80);

/* Set up time of next output, times of previous dump, log and periodic effect */
    if (*opflag == -1) {
	tout = *tstart;
    } else {
	n = (integer) ((d__1 = *time - *tstart, abs(d__1)) / *dtout) + 1;
	d__1 = (doublereal) n;
	d__2 = *tstop - *tstart;
	tout = *tstart + *dtout * d_sign(&d__1, &d__2);
	if ((*tstop - *tstart) * (tout - *tstop) > 0.) {
	    tout = *tstop;
	}
    }
    tdump = *time;
    tfun = *time;
    tlog = *time;

/* Convert to internal coordinates and velocities */
    (*coord)(time, &jcen[1], nbod, nbig, h0, &m[1], &xh[4], &vh[4], x, v, &
	    ngf[5], ngflag, &opt[1]);

/* ------------------------------------------------------------------------------ */

/*  MAIN  LOOP  STARTS  HERE */

L100:

/* Is it time for output ? */
    if ((d__1 = tout - *time, abs(d__1)) <= hby2 && *opflag >= -1) {

/* Beware: the integration may change direction at this point!!!! */
	if (*opflag == -1 && dtflag != 0) {
	    dtflag = 1;
	}

/* Convert to heliocentric coordinates and output data for all bodies */
	(*bcoord)(time, &jcen[1], nbod, nbig, h0, &m[1], x, v, &xh[4], &vh[4],
		 &ngf[5], ngflag, &opt[1]);
	mio_out__(time, &jcen[1], rcen, rmax, nbod, nbig, &m[1], &xh[4], &vh[
		4], &s[4], &rho[1], &stat[1], id + 8, &opt[1], opflag, algor, 
		outfile + 80, (ftnlen)8, (ftnlen)80);
	mio_ce__(time, tstart, rcen, rmax, nbod, nbig, &m[1], &stat[1], id + 
		8, &c__0, iclo, jclo, &opt[1], &stopflag, tclo, dclo, ixvclo, 
		jxvclo, mem + 80, &lmem[1], outfile + 80, &nstored, &c__0, (
		ftnlen)8, (ftnlen)80, (ftnlen)80);
	tmp0 = *tstop - tout;
/* Computing MIN */
	d__2 = abs(tmp0), d__3 = abs(*dtout);
	d__1 = min(d__2,d__3);
	tout += d_sign(&d__1, &tmp0);

/* Update the data dump files */
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    epoch[j - 1] = *time;
	}
	mio_dump__(time, tstart, tstop, dtout, algor, h0, tol, &jcen[1], rcen,
		 rmax, &en[1], &am[1], cefac, ndump, nfun, nbod, nbig, &m[1], 
		&xh[4], &vh[4], &s[4], &rho[1], &rceh[1], &stat[1], id + 8, &
		ngf[5], epoch, &opt[1], opflag, dumpfile + 80, mem + 80, &
		lmem[1], (ftnlen)8, (ftnlen)80, (ftnlen)80);
	tdump = *time;
    }

/* If integration has finished, convert to heliocentric coords and return */
    if ((d__1 = *tstop - *time, abs(d__1)) <= hby2 && *opflag >= 0) {
	(*bcoord)(time, &jcen[1], nbod, nbig, h0, &m[1], x, v, &xh[4], &vh[4],
		 &ngf[5], ngflag, &opt[1]);
	return 0;
    }

/* Make sure the integration is heading in the right direction */
L150:
    tmp0 = *tstop - *time;
    if (*opflag == -1) {
	tmp0 = *tstart - *time;
    }
    *h0 = d_sign(h0, &tmp0);

/* Save the current heliocentric coordinates and velocities */
    if (*algor == 1) {
	mco_iden__(time, &jcen[1], nbod, nbig, h0, &m[1], x, v, xh0, vh0, &
		ngf[5], ngflag, &opt[1]);
    } else {
	(*bcoord)(time, &jcen[1], nbod, nbig, h0, &m[1], x, v, xh0, vh0, &ngf[
		5], ngflag, &opt[1]);
    }
    (*onestep)(time, tstart, h0, tol, rmax, &en[1], &am[1], &jcen[1], rcen, 
	    nbod, nbig, &m[1], x, v, &s[4], rphys, rcrit, rce, &stat[1], id + 
	    8, &ngf[5], algor, &opt[1], &dtflag, ngflag, opflag, &colflag, &
	    nclo, iclo, jclo, dclo, tclo, ixvclo, jxvclo, outfile + 80, mem + 
	    80, &lmem[1], (ftnlen)8, (ftnlen)80, (ftnlen)80);
    *time += *h0;

/* ------------------------------------------------------------------------------ */

/*  CLOSE  ENCOUNTERS */

/* If encounter minima occurred, output details and decide whether to stop */
    if (nclo > 0 && *opflag >= -1) {
	itmp = 1;
	if (colflag != 0) {
	    itmp = 0;
	}
	mio_ce__(time, tstart, rcen, rmax, nbod, nbig, &m[1], &stat[1], id + 
		8, &nclo, iclo, jclo, &opt[1], &stopflag, tclo, dclo, ixvclo, 
		jxvclo, mem + 80, &lmem[1], outfile + 80, &nstored, &itmp, (
		ftnlen)8, (ftnlen)80, (ftnlen)80);
	if (stopflag == 1) {
	    return 0;
	}
    }

/* ------------------------------------------------------------------------------ */

/*  COLLISIONS */

/* If collisions occurred, output details and remove lost objects */
    if (colflag != 0) {

/* Reindex the surviving objects */
	(*bcoord)(time, &jcen[1], nbod, nbig, h0, &m[1], x, v, &xh[4], &vh[4],
		 &ngf[5], ngflag, &opt[1]);
	mxx_elim__(nbod, nbig, &m[1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1]
		, rcrit, &ngf[5], &stat[1], id + 8, mem + 80, &lmem[1], 
		outfile + 240, &itmp, (ftnlen)8, (ftnlen)80, (ftnlen)80);

/* Reset flags, and calculate new Hill radii and physical radii */
	dtflag = 1;
	if (*opflag >= 0) {
	    *opflag = 1;
	}
	mce_init__(tstart, algor, h0, &jcen[1], rcen, rmax, cefac, nbod, nbig,
		 &m[1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1], rphys, rce, 
		rcrit, id + 8, &opt[1], outfile + 160, &c__1, (ftnlen)8, (
		ftnlen)80);
	(*coord)(time, &jcen[1], nbod, nbig, h0, &m[1], &xh[4], &vh[4], x, v, 
		&ngf[5], ngflag, &opt[1]);
    }

/* ------------------------------------------------------------------------------ */

/*  COLLISIONS  WITH  CENTRAL  BODY */

/* Check for collisions with the central body */
    if (*algor == 1) {
	mco_iden__(time, &jcen[1], nbod, nbig, h0, &m[1], x, v, &xh[4], &vh[4]
		, &ngf[5], ngflag, &opt[1]);
    } else {
	(*bcoord)(time, &jcen[1], nbod, nbig, h0, &m[1], x, v, &xh[4], &vh[4],
		 &ngf[5], ngflag, &opt[1]);
    }
    itmp = 2;
    if (*algor == 11 || *algor == 12) {
	itmp = 3;
    }
    mce_cent__(time, h0, rcen, &jcen[1], &itmp, nbod, nbig, &m[1], xh0, vh0, &
	    xh[4], &vh[4], &nhit, jhit, thit, dhit, algor, &ngf[5], ngflag);

/* If something hit the central body, restore the coords prior to this step */
    if (nhit > 0) {
	mco_iden__(time, &jcen[1], nbod, nbig, h0, &m[1], xh0, vh0, &xh[4], &
		vh[4], &ngf[5], ngflag, &opt[1]);
	*time -= *h0;

/* Merge the object(s) with the central body */
	i__1 = nhit;
	for (k = 1; k <= i__1; ++k) {
	    i__ = 1;
	    j = jhit[k - 1];
	    mce_coll__(&thit[k - 1], tstart, &en[3], &jcen[1], &i__, &j, nbod,
		     nbig, &m[1], &xh[4], &vh[4], &s[4], rphys, &stat[1], id 
		    + 8, &opt[1], mem + 80, &lmem[1], outfile + 240, (ftnlen)
		    8, (ftnlen)80, (ftnlen)80);
	}

/* Remove lost objects, reset flags and recompute Hill and physical radii */
	mxx_elim__(nbod, nbig, &m[1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1]
		, rcrit, &ngf[5], &stat[1], id + 8, mem + 80, &lmem[1], 
		outfile + 240, &itmp, (ftnlen)8, (ftnlen)80, (ftnlen)80);
	if (*opflag >= 0) {
	    *opflag = 1;
	}
	dtflag = 1;
	mce_init__(tstart, algor, h0, &jcen[1], rcen, rmax, cefac, nbod, nbig,
		 &m[1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1], rphys, rce, 
		rcrit, id + 8, &opt[1], outfile + 160, &c__0, (ftnlen)8, (
		ftnlen)80);
	if (*algor == 1) {
	    mco_iden__(time, &jcen[1], nbod, nbig, h0, &m[1], &xh[4], &vh[4], 
		    x, v, &ngf[5], ngflag, &opt[1]);
	} else {
	    (*coord)(time, &jcen[1], nbod, nbig, h0, &m[1], &xh[4], &vh[4], x,
		     v, &ngf[5], ngflag, &opt[1]);
	}

/* Redo that integration time step */
	goto L150;
    }

/* ------------------------------------------------------------------------------ */

/*  DATA  DUMP  AND  PROGRESS  REPORT */

/* Convert to heliocentric coords and do the data dump */
    if ((d__1 = *time - tdump, abs(d__1)) >= abs(dtdump) && *opflag >= -1) {
	(*bcoord)(time, &jcen[1], nbod, nbig, h0, &m[1], x, v, &xh[4], &vh[4],
		 &ngf[5], ngflag, &opt[1]);
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    epoch[j - 1] = *time;
	}
	mio_ce__(time, tstart, rcen, rmax, nbod, nbig, &m[1], &stat[1], id + 
		8, &c__0, iclo, jclo, &opt[1], &stopflag, tclo, dclo, ixvclo, 
		jxvclo, mem + 80, &lmem[1], outfile + 80, &nstored, &c__0, (
		ftnlen)8, (ftnlen)80, (ftnlen)80);
	mio_dump__(time, tstart, tstop, dtout, algor, h0, tol, &jcen[1], rcen,
		 rmax, &en[1], &am[1], cefac, ndump, nfun, nbod, nbig, &m[1], 
		&xh[4], &vh[4], &s[4], &rho[1], &rceh[1], &stat[1], id + 8, &
		ngf[5], epoch, &opt[1], opflag, dumpfile + 80, mem + 80, &
		lmem[1], (ftnlen)8, (ftnlen)80, (ftnlen)80);
	tdump = *time;
    }

/* Convert to heliocentric coords and write a progress report to the log file */
    if ((d__1 = *time - tlog, abs(d__1)) >= abs(dtdump) && *opflag >= 0) {
	(*bcoord)(time, &jcen[1], nbod, nbig, h0, &m[1], x, v, &xh[4], &vh[4],
		 &ngf[5], ngflag, &opt[1]);
	mxx_en__(&jcen[1], nbod, nbig, &m[1], &xh[4], &vh[4], &s[4], &en[2], &
		am[2]);
	mio_log__(time, tstart, &en[1], &am[1], &opt[1], mem + 80, &lmem[1], (
		ftnlen)80);
	tlog = *time;
    }

/* ------------------------------------------------------------------------------ */

/*  CHECK  FOR  EJECTIONS  AND  DO  OTHER  PERIODIC  EFFECTS */

    if ((d__1 = *time - tfun, abs(d__1)) >= abs(dtfun) && *opflag >= -1) {
	if (*algor == 1) {
	    mco_iden__(time, &jcen[1], nbod, nbig, h0, &m[1], x, v, &xh[4], &
		    vh[4], &ngf[5], ngflag, &opt[1]);
	} else {
	    (*bcoord)(time, &jcen[1], nbod, nbig, h0, &m[1], x, v, &xh[4], &
		    vh[4], &ngf[5], ngflag, &opt[1]);
	}

/* Recompute close encounter limits, to allow for changes in Hill radii */
	mce_hill__(nbod, &m[1], &xh[4], &vh[4], rce, a);
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    rce[j - 1] *= rceh[j];
	}

/* Check for ejections */
	itmp = 2;
	if (*algor == 11 || *algor == 12) {
	    itmp = 3;
	}
	mxx_ejec__(time, tstart, rmax, &en[1], &am[1], &jcen[1], &itmp, nbod, 
		nbig, &m[1], &xh[4], &vh[4], &s[4], &stat[1], id + 8, &opt[1],
		 &ejflag, outfile + 240, mem + 80, &lmem[1], (ftnlen)8, (
		ftnlen)80, (ftnlen)80);

/* Remove ejected objects, reset flags, calculate new Hill and physical radii */
	if (ejflag != 0) {
	    mxx_elim__(nbod, nbig, &m[1], &xh[4], &vh[4], &s[4], &rho[1], &
		    rceh[1], rcrit, &ngf[5], &stat[1], id + 8, mem + 80, &
		    lmem[1], outfile + 240, &itmp, (ftnlen)8, (ftnlen)80, (
		    ftnlen)80);
	    if (*opflag >= 0) {
		*opflag = 1;
	    }
	    dtflag = 1;
	    mce_init__(tstart, algor, h0, &jcen[1], rcen, rmax, cefac, nbod, 
		    nbig, &m[1], &xh[4], &vh[4], &s[4], &rho[1], &rceh[1], 
		    rphys, rce, rcrit, id + 8, &opt[1], outfile + 160, &c__0, 
		    (ftnlen)8, (ftnlen)80);
	    if (*algor == 1) {
		mco_iden__(time, &jcen[1], nbod, nbig, h0, &m[1], &xh[4], &vh[
			4], x, v, &ngf[5], ngflag, &opt[1]);
	    } else {
		(*coord)(time, &jcen[1], nbod, nbig, h0, &m[1], &xh[4], &vh[4]
			, x, v, &ngf[5], ngflag, &opt[1]);
	    }
	}
	tfun = *time;
    }

/* Go on to the next time step */
    goto L100;

/* ------------------------------------------------------------------------------ */

} /* mal_hcon__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*    MCE_BOX.FOR    (ErikSoft   30 September 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Given initial and final coordinates and velocities, the routine returns */
/* the X and Y coordinates of a box bounding the motion in between the */
/* end points. */

/* If the X or Y velocity changes sign, the routine performs a quadratic */
/* interpolation to estimate the corresponding extreme value of X or Y. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mce_box__(integer *nbod, doublereal *h__, doublereal *x0,
	 doublereal *v0, doublereal *x1, doublereal *v1, doublereal *xmin, 
	doublereal *xmax, doublereal *ymin, doublereal *ymax)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j;
    static doublereal temp;



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
    --ymax;
    --ymin;
    --xmax;
    --xmin;
    v1 -= 4;
    x1 -= 4;
    v0 -= 4;
    x0 -= 4;

    /* Function Body */
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
	d__1 = x0[j * 3 + 1], d__2 = x1[j * 3 + 1];
	xmin[j] = min(d__1,d__2);
/* Computing MAX */
	d__1 = x0[j * 3 + 1], d__2 = x1[j * 3 + 1];
	xmax[j] = max(d__1,d__2);
/* Computing MIN */
	d__1 = x0[j * 3 + 2], d__2 = x1[j * 3 + 2];
	ymin[j] = min(d__1,d__2);
/* Computing MAX */
	d__1 = x0[j * 3 + 2], d__2 = x1[j * 3 + 2];
	ymax[j] = max(d__1,d__2);

/* If velocity changes sign, do an interpolation */
	if (v0[j * 3 + 1] < 0. && v1[j * 3 + 1] > 0. || v0[j * 3 + 1] > 0. && 
		v1[j * 3 + 1] < 0.) {
	    temp = (v0[j * 3 + 1] * x1[j * 3 + 1] - v1[j * 3 + 1] * x0[j * 3 
		    + 1] - *h__ * .5 * v0[j * 3 + 1] * v1[j * 3 + 1]) / (v0[j 
		    * 3 + 1] - v1[j * 3 + 1]);
/* Computing MIN */
	    d__1 = xmin[j];
	    xmin[j] = min(d__1,temp);
/* Computing MAX */
	    d__1 = xmax[j];
	    xmax[j] = max(d__1,temp);
	}

	if (v0[j * 3 + 2] < 0. && v1[j * 3 + 2] > 0. || v0[j * 3 + 2] > 0. && 
		v1[j * 3 + 2] < 0.) {
	    temp = (v0[j * 3 + 2] * x1[j * 3 + 2] - v1[j * 3 + 2] * x0[j * 3 
		    + 2] - *h__ * .5 * v0[j * 3 + 2] * v1[j * 3 + 2]) / (v0[j 
		    * 3 + 2] - v1[j * 3 + 2]);
/* Computing MIN */
	    d__1 = ymin[j];
	    ymin[j] = min(d__1,temp);
/* Computing MAX */
	    d__1 = ymax[j];
	    ymax[j] = max(d__1,temp);
	}
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mce_box__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*    MCE_CENT.FOR    (ErikSoft   4 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Checks all objects with index I >= I0, to see if they have had a collision */
/* with the central body in a time interval H, when given the initial and */
/* final coordinates and velocities. The routine uses cubic interpolation */
/* to estimate the minimum separations. */

/* N.B. All coordinates & velocities must be with respect to the central body!! */
/* === */
/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mce_cent__(doublereal *time, doublereal *h__, doublereal 
	*rcen, doublereal *jcen, integer *i0, integer *nbod, integer *nbig, 
	doublereal *m, doublereal *x0, doublereal *v0, doublereal *x1, 
	doublereal *v1, integer *nhit, integer *jhit, doublereal *thit, 
	doublereal *dhit, integer *algor, doublereal *ngf, integer *ngflag)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal), d_sign(doublereal *, 
	    doublereal *), sin(doublereal), d_mod(doublereal *, doublereal *),
	     sinh(doublereal);

    /* Local variables */
    extern doublereal mco_acsh__(doublereal *);
    static doublereal a, e;
    static integer j;
    static doublereal p, q, h2, m0, r0, u0, v2, mm, hx, hy, hz, rr0, rr1, rv0,
	     rv1, vu0[6000]	/* was [3][2000] */, vu1[6000]	/* was [3][
	    2000] */, xu0[6000]	/* was [3][2000] */, xu1[6000]	/* was [3][
	    2000] */, mcen, mhit, temp, uhit, rcen2;



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
    --jcen;
    ngf -= 5;
    v1 -= 4;
    x1 -= 4;
    v0 -= 4;
    x0 -= 4;
    --m;
    --jhit;
    --thit;
    --dhit;

    /* Function Body */
    if (*i0 <= 0) {
	*i0 = 2;
    }
    *nhit = 0;
    rcen2 = *rcen * *rcen;
    mcen = m[1];

/* If using close-binary code, convert to coords with respect to the binary */
/*      if (algor.eq.11) then */
/*        mcen = m(1) + m(2) */
/*        call mco_h2ub (temp,jcen,nbod,nbig,h,m,x0,v0,xu0,vu0,ngf,ngflag) */
/*        call mco_h2ub (temp,jcen,nbod,nbig,h,m,x1,v1,xu1,vu1,ngf,ngflag) */
/*      end if */

/* Check for collisions with the central body */
    i__1 = *nbod;
    for (j = *i0; j <= i__1; ++j) {
	if (*algor == 11) {
	    rr0 = xu0[j * 3 - 3] * xu0[j * 3 - 3] + xu0[j * 3 - 2] * xu0[j * 
		    3 - 2] + xu0[j * 3 - 1] * xu0[j * 3 - 1];
	    rr1 = xu1[j * 3 - 3] * xu1[j * 3 - 3] + xu1[j * 3 - 2] * xu1[j * 
		    3 - 2] + xu1[j * 3 - 1] * xu1[j * 3 - 1];
	    rv0 = vu0[j * 3 - 3] * xu0[j * 3 - 3] + vu0[j * 3 - 2] * xu0[j * 
		    3 - 2] + vu0[j * 3 - 1] * xu0[j * 3 - 1];
	    rv1 = vu1[j * 3 - 3] * xu1[j * 3 - 3] + vu1[j * 3 - 2] * xu1[j * 
		    3 - 2] + vu1[j * 3 - 1] * xu1[j * 3 - 1];
	} else {
	    rr0 = x0[j * 3 + 1] * x0[j * 3 + 1] + x0[j * 3 + 2] * x0[j * 3 + 
		    2] + x0[j * 3 + 3] * x0[j * 3 + 3];
	    rr1 = x1[j * 3 + 1] * x1[j * 3 + 1] + x1[j * 3 + 2] * x1[j * 3 + 
		    2] + x1[j * 3 + 3] * x1[j * 3 + 3];
	    rv0 = v0[j * 3 + 1] * x0[j * 3 + 1] + v0[j * 3 + 2] * x0[j * 3 + 
		    2] + v0[j * 3 + 3] * x0[j * 3 + 3];
	    rv1 = v1[j * 3 + 1] * x1[j * 3 + 1] + v1[j * 3 + 2] * x1[j * 3 + 
		    2] + v1[j * 3 + 3] * x1[j * 3 + 3];
	}

/* If inside the central body, or passing through pericentre, use 2-body approx. */
	if (rv0 * *h__ <= 0. && rv1 * *h__ >= 0. || min(rr0,rr1) <= rcen2) {
	    if (*algor == 11) {
		hx = xu0[j * 3 - 2] * vu0[j * 3 - 1] - xu0[j * 3 - 1] * vu0[j 
			* 3 - 2];
		hy = xu0[j * 3 - 1] * vu0[j * 3 - 3] - xu0[j * 3 - 3] * vu0[j 
			* 3 - 1];
		hz = xu0[j * 3 - 3] * vu0[j * 3 - 2] - xu0[j * 3 - 2] * vu0[j 
			* 3 - 3];
		v2 = vu0[j * 3 - 3] * vu0[j * 3 - 3] + vu0[j * 3 - 2] * vu0[j 
			* 3 - 2] + vu0[j * 3 - 1] * vu0[j * 3 - 1];
	    } else {
		hx = x0[j * 3 + 2] * v0[j * 3 + 3] - x0[j * 3 + 3] * v0[j * 3 
			+ 2];
		hy = x0[j * 3 + 3] * v0[j * 3 + 1] - x0[j * 3 + 1] * v0[j * 3 
			+ 3];
		hz = x0[j * 3 + 1] * v0[j * 3 + 2] - x0[j * 3 + 2] * v0[j * 3 
			+ 1];
		v2 = v0[j * 3 + 1] * v0[j * 3 + 1] + v0[j * 3 + 2] * v0[j * 3 
			+ 2] + v0[j * 3 + 3] * v0[j * 3 + 3];
	    }
	    h2 = hx * hx + hy * hy + hz * hz;
	    p = h2 / (mcen + m[j]);
	    r0 = sqrt(rr0);
	    temp = p * (v2 / (mcen + m[j]) - 2. / r0) + 1.;
	    e = sqrt((max(temp,0.)));
	    q = p / (e + 1.);

/* If the object hit the central body */
	    if (q <= *rcen) {
		++(*nhit);
		jhit[*nhit] = j;
		dhit[*nhit] = *rcen;

/* Time of impact relative to the end of the timestep */
		if (e < 1.) {
		    a = q / (1. - e);
		    d__1 = acos((1. - *rcen / a) / e);
		    d__2 = -(*h__);
		    uhit = d_sign(&d__1, &d__2);
		    d__1 = acos((1. - r0 / a) / e);
		    u0 = d_sign(&d__1, &rv0);
		    d__1 = uhit - e * sin(uhit) + 3.141592653589793;
		    mhit = d_mod(&d__1, &c_b58) - 3.141592653589793;
		    d__1 = u0 - e * sin(u0) + 3.141592653589793;
		    m0 = d_mod(&d__1, &c_b58) - 3.141592653589793;
		} else {
		    a = q / (e - 1.);
		    d__2 = (1. - *rcen / a) / e;
		    d__1 = mco_acsh__(&d__2);
		    d__3 = -(*h__);
		    uhit = d_sign(&d__1, &d__3);
		    d__2 = (1. - r0 / a) / e;
		    d__1 = mco_acsh__(&d__2);
		    u0 = d_sign(&d__1, &rv0);
		    d__1 = uhit - e * sinh(uhit) + 3.141592653589793;
		    mhit = d_mod(&d__1, &c_b58) - 3.141592653589793;
		    d__1 = u0 - e * sinh(u0) + 3.141592653589793;
		    m0 = d_mod(&d__1, &c_b58) - 3.141592653589793;
		}
		mm = sqrt((mcen + m[j]) / (a * a * a));
		thit[*nhit] = (mhit - m0) / mm + *time;
	    }
	}
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mce_cent__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCE_COLL.FOR    (ErikSoft   2 October 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Resolves a collision between two objects, using the collision model chosen */
/* by the user. Also writes a message to the information file, and updates the */
/* value of ELOST, the change in energy due to collisions and ejections. */

/* N.B. All coordinates and velocities must be with respect to central body. */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mce_coll__(doublereal *time, doublereal *tstart, 
	doublereal *elost, doublereal *jcen, integer *i__, integer *j, 
	integer *nbod, integer *nbig, doublereal *m, doublereal *xh, 
	doublereal *vh, doublereal *s, doublereal *rphys, integer *stat, char 
	*id, integer *opt, char *mem, integer *lmem, char *outfile, ftnlen 
	id_len, ftnlen mem_len, ftnlen outfile_len)
{
    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int mce_merg__(doublereal *, integer *, integer *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *);
    static doublereal t1;
    static char fcol[38];
    static integer year, itmp, month;
    static char flost[38], tstring[6];
    extern /* Subroutine */ int mio_jd2y__(doublereal *, integer *, integer *,
	     doublereal *);

    /* Fortran I/O blocks */
    static cilist io___164 = { 0, 23, 0, flost, 0 };
    static cilist io___166 = { 0, 23, 0, fcol, 0 };
    static cilist io___168 = { 0, 23, 0, flost, 0 };
    static cilist io___169 = { 0, 23, 0, fcol, 0 };




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

/* If two bodies collided, check that the less massive one is removed */
/* (unless the more massive one is a Small body) */
    /* Parameter adjustments */
    --jcen;
    id -= 8;
    --stat;
    --rphys;
    s -= 4;
    vh -= 4;
    xh -= 4;
    --m;
    --opt;
    mem -= 80;
    --lmem;

    /* Function Body */
    if (*i__ != 0) {
	if (m[*j] > m[*i__] && *j <= *nbig) {
	    itmp = *i__;
	    *i__ = *j;
	    *j = itmp;
	}
    }

/* Write message to info file (I=0 implies collision with the central body) */
L10:
    o__1.oerr = 1;
    o__1.ounit = 23;
    o__1.ofnmlen = 80;
    o__1.ofnm = outfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = "append";
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L10;
    }

    if (opt[3] == 1) {
	mio_jd2y__(time, &year, &month, &t1);
	if (*i__ == 0) {
	    s_copy(flost, "(1x,a8,a,i10,1x,i2,1x,f8.5)", (ftnlen)38, (ftnlen)
		    27);
	    s_wsfe(&io___164);
	    do_fio(&c__1, id + (*j << 3), (ftnlen)8);
	    do_fio(&c__1, mem + 5360, lmem[67]);
	    do_fio(&c__1, (char *)&year, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&month, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    s_copy(fcol, "(1x,a8,a,a8,a,i10,1x,i2,1x,f4.1)", (ftnlen)38, (
		    ftnlen)32);
	    s_wsfe(&io___166);
	    do_fio(&c__1, id + (*i__ << 3), (ftnlen)8);
	    do_fio(&c__1, mem + 5520, lmem[69]);
	    do_fio(&c__1, id + (*j << 3), (ftnlen)8);
	    do_fio(&c__1, mem + 5680, lmem[71]);
	    do_fio(&c__1, (char *)&year, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&month, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    } else {
	if (opt[3] == 3) {
	    t1 = (*time - *tstart) / 365.25;
	    s_copy(tstring, mem + 160, (ftnlen)6, (ftnlen)80);
	    s_copy(flost, "(1x,a8,a,f18.7,a)", (ftnlen)38, (ftnlen)17);
	    s_copy(fcol, "(1x,a8,a,a8,a,1x,f14.3,a)", (ftnlen)38, (ftnlen)25);
	} else {
	    if (opt[3] == 0) {
		t1 = *time;
	    }
	    if (opt[3] == 2) {
		t1 = *time - *tstart;
	    }
	    s_copy(tstring, mem + 80, (ftnlen)6, (ftnlen)80);
	    s_copy(flost, "(1x,a8,a,f18.5,a)", (ftnlen)38, (ftnlen)17);
	    s_copy(fcol, "(1x,a8,a,a8,a,1x,f14.1,a)", (ftnlen)38, (ftnlen)25);
	}
	if (*i__ == 0 || *i__ == 1) {
	    s_wsfe(&io___168);
	    do_fio(&c__1, id + (*j << 3), (ftnlen)8);
	    do_fio(&c__1, mem + 5360, lmem[67]);
	    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, tstring, (ftnlen)6);
	    e_wsfe();
	} else {
	    s_wsfe(&io___169);
	    do_fio(&c__1, id + (*i__ << 3), (ftnlen)8);
	    do_fio(&c__1, mem + 5520, lmem[69]);
	    do_fio(&c__1, id + (*j << 3), (ftnlen)8);
	    do_fio(&c__1, mem + 5680, lmem[71]);
	    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, tstring, (ftnlen)6);
	    e_wsfe();
	}
    }
    cl__1.cerr = 0;
    cl__1.cunit = 23;
    cl__1.csta = 0;
    f_clos(&cl__1);

/* Do the collision (inelastic merger) */
    mce_merg__(&jcen[1], i__, j, nbod, nbig, &m[1], &xh[4], &vh[4], &s[4], &
	    stat[1], elost);

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mce_coll__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCE_HILL.FOR    (ErikSoft   4 October 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates the Hill radii for all objects given their masses, M, */
/* coordinates, X, and velocities, V; plus the mass of the central body, M(1) */
/* Where HILL = a * (m/3*m(1))^(1/3) */

/* If the orbit is hyperbolic or parabolic, the Hill radius is calculated using: */
/*       HILL = r * (m/3*m(1))^(1/3) */
/* where R is the current distance from the central body. */

/* The routine also gives the semi-major axis, A, of each object's orbit. */

/* N.B. Designed to use heliocentric coordinates, but should be adequate using */
/* ===  barycentric coordinates. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mce_hill__(integer *nbod, doublereal *m, doublereal *x, 
	doublereal *v, doublereal *hill, doublereal *a)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer j;
    static doublereal r__, v2, gm;
    extern /* Subroutine */ int mco_x2a__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


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



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    --a;
    --hill;
    v -= 4;
    x -= 4;
    --m;

    /* Function Body */
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	gm = m[1] + m[j];
	mco_x2a__(&gm, &x[j * 3 + 1], &x[j * 3 + 2], &x[j * 3 + 3], &v[j * 3 
		+ 1], &v[j * 3 + 2], &v[j * 3 + 3], &a[j], &r__, &v2);
/* If orbit is hyperbolic, use the distance rather than the semi-major axis */
	if (a[j] <= 0.) {
	    a[j] = r__;
	}
	d__1 = m[j] * .3333333333333333 / m[1];
	hill[j] = a[j] * pow_dd(&d__1, &c_b95);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mce_hill__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCE_INIT.FOR    (ErikSoft   28 February 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates close-approach limits RCE (in AU) and physical radii RPHYS */
/* (in AU) for all objects, given their masses M, coordinates X, velocities */
/* V, densities RHO, and close-approach limits RCEH (in Hill radii). */

/* Also calculates the changeover distance RCRIT, used by the hybrid */
/* symplectic integrator. RCRIT is defined to be the larger of N1*HILL and */
/* N2*H*VMAX, where HILL is the Hill radius, H is the timestep, VMAX is the */
/* largest expected velocity of any body, and N1, N2 are parameters (see */
/* section 4.2 of Chambers 1999, Monthly Notices, vol 304, p793-799). */

/* N.B. Designed to use heliocentric coordinates, but should be adequate using */
/* ===  barycentric coordinates. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mce_init__(doublereal *tstart, integer *algor, 
	doublereal *h__, doublereal *jcen, doublereal *rcen, doublereal *rmax,
	 doublereal *cefac, integer *nbod, integer *nbig, doublereal *m, 
	doublereal *x, doublereal *v, doublereal *s, doublereal *rho, 
	doublereal *rceh, doublereal *rphys, doublereal *rce, doublereal *
	rcrit, char *id, integer *opt, char *outfile, integer *rcritflag, 
	ftnlen id_len, ftnlen outfile_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    char ch__1[8];
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int mce_hill__(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal a[2000];
    static char c__[80*2000];
    static integer j;
    static doublereal k_2__, amin, hill[2000], temp, vmax, rcen_2__;
    static char header[80];
    static doublereal rhocgs;
    extern /* Character */ VOID mio_fl2c__(char *, ftnlen, doublereal *), 
	    mio_re2c__(char *, ftnlen, doublereal *, doublereal *, doublereal 
	    *);

    /* Fortran I/O blocks */
    static cilist io___185 = { 0, 22, 0, "(a1,a2,i2,a62,i1)", 0 };
    static cilist io___186 = { 0, 22, 0, "(a51)", 0 };




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



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    --jcen;
    id -= 8;
    --rcrit;
    --rce;
    --rphys;
    --rceh;
    --rho;
    s -= 4;
    v -= 4;
    x -= 4;
    --m;
    --opt;

    /* Function Body */
    rhocgs = 498.06095345055081;
    k_2__ = 3379.3806811609443;
    rcen_2__ = 1. / (*rcen * *rcen);
    amin = 9.9e29;

/* Calculate the Hill radii */
    mce_hill__(nbod, &m[1], &x[4], &v[4], hill, a);

/* Determine the maximum close-encounter distances, and the physical radii */
    temp = m[1] * 2.25 / 3.141592653589793;
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	rce[j] = hill[j - 1] * rceh[j];
	d__1 = temp / rho[j];
	rphys[j] = hill[j - 1] / a[j - 1] * pow_dd(&d__1, &c_b95);
/* Computing MIN */
	d__1 = a[j - 1];
	amin = min(d__1,amin);
    }

/* If required, calculate the changeover distance used by hybrid algorithm */
    if (*rcritflag == 1) {
	vmax = sqrt(m[1] / amin);
	temp = *h__ * .4 * vmax;
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
/* Computing MAX */
	    d__1 = hill[j - 1] * *cefac;
	    rcrit[j] = max(d__1,temp);
	}
    }

/* Write list of object's identities to close-encounter output file */
    mio_fl2c__(ch__1, (ftnlen)8, tstart);
    s_copy(header, ch__1, (ftnlen)8, (ftnlen)8);
    d__1 = (doublereal) (*nbig - 1);
    mio_re2c__(ch__1, (ftnlen)8, &d__1, &c_b98, &c_b99);
    s_copy(header + 8, ch__1, (ftnlen)8, (ftnlen)8);
    d__1 = (doublereal) (*nbod - *nbig);
    mio_re2c__(ch__1, (ftnlen)8, &d__1, &c_b98, &c_b99);
    s_copy(header + 11, ch__1, (ftnlen)8, (ftnlen)8);
    d__1 = m[1] * k_2__;
    mio_fl2c__(ch__1, (ftnlen)8, &d__1);
    s_copy(header + 14, ch__1, (ftnlen)8, (ftnlen)8);
    d__1 = jcen[1] * rcen_2__;
    mio_fl2c__(ch__1, (ftnlen)8, &d__1);
    s_copy(header + 22, ch__1, (ftnlen)8, (ftnlen)8);
    d__1 = jcen[2] * rcen_2__ * rcen_2__;
    mio_fl2c__(ch__1, (ftnlen)8, &d__1);
    s_copy(header + 30, ch__1, (ftnlen)8, (ftnlen)8);
    d__1 = jcen[3] * rcen_2__ * rcen_2__ * rcen_2__;
    mio_fl2c__(ch__1, (ftnlen)8, &d__1);
    s_copy(header + 38, ch__1, (ftnlen)8, (ftnlen)8);
    mio_fl2c__(ch__1, (ftnlen)8, rcen);
    s_copy(header + 46, ch__1, (ftnlen)8, (ftnlen)8);
    mio_fl2c__(ch__1, (ftnlen)8, rmax);
    s_copy(header + 54, ch__1, (ftnlen)8, (ftnlen)8);

    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	d__1 = (doublereal) (j - 1);
	mio_re2c__(ch__1, (ftnlen)8, &d__1, &c_b98, &c_b99);
	s_copy(c__ + (j - 1) * 80, ch__1, (ftnlen)8, (ftnlen)8);
	s_copy(c__ + ((j - 1) * 80 + 3), id + (j << 3), (ftnlen)8, (ftnlen)8);
	d__1 = m[j] * k_2__;
	mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	s_copy(c__ + ((j - 1) * 80 + 11), ch__1, (ftnlen)8, (ftnlen)8);
	d__1 = s[j * 3 + 1] * k_2__;
	mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	s_copy(c__ + ((j - 1) * 80 + 19), ch__1, (ftnlen)8, (ftnlen)8);
	d__1 = s[j * 3 + 2] * k_2__;
	mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	s_copy(c__ + ((j - 1) * 80 + 27), ch__1, (ftnlen)8, (ftnlen)8);
	d__1 = s[j * 3 + 3] * k_2__;
	mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	s_copy(c__ + ((j - 1) * 80 + 35), ch__1, (ftnlen)8, (ftnlen)8);
	d__1 = rho[j] / rhocgs;
	mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	s_copy(c__ + ((j - 1) * 80 + 43), ch__1, (ftnlen)8, (ftnlen)8);
    }

/* Write compressed output to file */
L50:
    o__1.oerr = 1;
    o__1.ounit = 22;
    o__1.ofnmlen = 80;
    o__1.ofnm = outfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = "append";
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L50;
    }
    s_wsfe(&io___185);
    do_fio(&c__1, "\f", (ftnlen)1);
    do_fio(&c__1, "6a", (ftnlen)2);
    do_fio(&c__1, (char *)&(*algor), (ftnlen)sizeof(integer));
    do_fio(&c__1, header, (ftnlen)62);
    do_fio(&c__1, (char *)&opt[4], (ftnlen)sizeof(integer));
    e_wsfe();
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	s_wsfe(&io___186);
	do_fio(&c__1, c__ + (j - 1) * 80, (ftnlen)51);
	e_wsfe();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 22;
    cl__1.csta = 0;
    f_clos(&cl__1);

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mce_init__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCE_MERG.FOR    (ErikSoft   2 October 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% c */
/* Author: John E. Chambers */

/* Merges objects I and J inelastically to produce a single new body by */
/* conserving mass and linear momentum. */
/*   If J <= NBIG, then J is a Big body */
/*   If J >  NBIG, then J is a Small body */
/*   If I = 0, then I is the central body */

/* N.B. All coordinates and velocities must be with respect to central body. */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mce_merg__(doublereal *jcen, integer *i__, integer *j, 
	integer *nbod, integer *nbig, doublereal *m, doublereal *xh, 
	doublereal *vh, doublereal *s, integer *stat, doublereal *elost)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k;
    static doublereal e0, e1, l2, du, dv, dx, dy, dz, dw, tmp1, tmp2, msum, 
	    mredu, msum_1__;
    extern /* Subroutine */ int mxx_en__(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
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

/* If a body hits the central body */
    /* Parameter adjustments */
    --jcen;
    --stat;
    s -= 4;
    vh -= 4;
    xh -= 4;
    --m;

    /* Function Body */
    if (*i__ <= 1) {
	mxx_en__(&jcen[1], nbod, nbig, &m[1], &xh[4], &vh[4], &s[4], &e0, &l2)
		;

/* If a body hit the central body... */
	msum = m[1] + m[*j];
	msum_1__ = 1. / msum;
	mredu = m[1] * m[*j] * msum_1__;
	dx = xh[*j * 3 + 1];
	dy = xh[*j * 3 + 2];
	dz = xh[*j * 3 + 3];
	du = vh[*j * 3 + 1];
	dv = vh[*j * 3 + 2];
	dw = vh[*j * 3 + 3];

/* Calculate new spin angular momentum of the central body */
	s[4] = s[4] + s[*j * 3 + 1] + mredu * (dy * dw - dz * dv);
	s[5] = s[5] + s[*j * 3 + 2] + mredu * (dz * du - dx * dw);
	s[6] = s[6] + s[*j * 3 + 3] + mredu * (dx * dv - dy * du);

/* Calculate shift in barycentric coords and velocities of central body */
	tmp2 = m[*j] * msum_1__;
	xh[4] = tmp2 * xh[*j * 3 + 1];
	xh[5] = tmp2 * xh[*j * 3 + 2];
	xh[6] = tmp2 * xh[*j * 3 + 3];
	vh[4] = tmp2 * vh[*j * 3 + 1];
	vh[5] = tmp2 * vh[*j * 3 + 2];
	vh[6] = tmp2 * vh[*j * 3 + 3];
	m[1] = msum;
	m[*j] = 0.;
	s[*j * 3 + 1] = 0.;
	s[*j * 3 + 2] = 0.;
	s[*j * 3 + 3] = 0.;

/* Shift the heliocentric coordinates and velocities of all bodies */
	i__1 = *nbod;
	for (k = 2; k <= i__1; ++k) {
	    xh[k * 3 + 1] -= xh[4];
	    xh[k * 3 + 2] -= xh[5];
	    xh[k * 3 + 3] -= xh[6];
	    vh[k * 3 + 1] -= vh[4];
	    vh[k * 3 + 2] -= vh[5];
	    vh[k * 3 + 3] -= vh[6];
	}

/* Calculate energy loss due to the collision */
	mxx_en__(&jcen[1], nbod, nbig, &m[1], &xh[4], &vh[4], &s[4], &e1, &l2)
		;
	*elost += e0 - e1;
    } else {

/* If two bodies collided... */
	msum = m[*i__] + m[*j];
	msum_1__ = 1. / msum;
	mredu = m[*i__] * m[*j] * msum_1__;
	dx = xh[*i__ * 3 + 1] - xh[*j * 3 + 1];
	dy = xh[*i__ * 3 + 2] - xh[*j * 3 + 2];
	dz = xh[*i__ * 3 + 3] - xh[*j * 3 + 3];
	du = vh[*i__ * 3 + 1] - vh[*j * 3 + 1];
	dv = vh[*i__ * 3 + 2] - vh[*j * 3 + 2];
	dw = vh[*i__ * 3 + 3] - vh[*j * 3 + 3];

/* Calculate energy loss due to the collision */
	*elost = *elost + mredu * .5 * (du * du + dv * dv + dw * dw) - m[*i__]
		 * m[*j] / sqrt(dx * dx + dy * dy + dz * dz);

/* Calculate spin angular momentum of the new body */
	s[*i__ * 3 + 1] = s[*i__ * 3 + 1] + s[*j * 3 + 1] + mredu * (dy * dw 
		- dz * dv);
	s[*i__ * 3 + 2] = s[*i__ * 3 + 2] + s[*j * 3 + 2] + mredu * (dz * du 
		- dx * dw);
	s[*i__ * 3 + 3] = s[*i__ * 3 + 3] + s[*j * 3 + 3] + mredu * (dx * dv 
		- dy * du);

/* Calculate new coords and velocities by conserving centre of mass & momentum */
	tmp1 = m[*i__] * msum_1__;
	tmp2 = m[*j] * msum_1__;
	xh[*i__ * 3 + 1] = xh[*i__ * 3 + 1] * tmp1 + xh[*j * 3 + 1] * tmp2;
	xh[*i__ * 3 + 2] = xh[*i__ * 3 + 2] * tmp1 + xh[*j * 3 + 2] * tmp2;
	xh[*i__ * 3 + 3] = xh[*i__ * 3 + 3] * tmp1 + xh[*j * 3 + 3] * tmp2;
	vh[*i__ * 3 + 1] = vh[*i__ * 3 + 1] * tmp1 + vh[*j * 3 + 1] * tmp2;
	vh[*i__ * 3 + 2] = vh[*i__ * 3 + 2] * tmp1 + vh[*j * 3 + 2] * tmp2;
	vh[*i__ * 3 + 3] = vh[*i__ * 3 + 3] * tmp1 + vh[*j * 3 + 3] * tmp2;
	m[*i__] = msum;
    }

/* Flag the lost body for removal, and move it away from the new body */
    stat[*j] = -2;
    xh[*j * 3 + 1] = -xh[*j * 3 + 1];
    xh[*j * 3 + 2] = -xh[*j * 3 + 2];
    xh[*j * 3 + 3] = -xh[*j * 3 + 3];
    vh[*j * 3 + 1] = -vh[*j * 3 + 1];
    vh[*j * 3 + 2] = -vh[*j * 3 + 2];
    vh[*j * 3 + 3] = -vh[*j * 3 + 3];
    m[*j] = 0.;
    s[*j * 3 + 1] = 0.;
    s[*j * 3 + 2] = 0.;
    s[*j * 3 + 3] = 0.;

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mce_merg__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCE_MIN.FOR    (ErikSoft  1 December 1998) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates minimum value of a quantity D, within an interval H, given initial */
/* and final values D0, D1, and their derivatives D0T, D1T, using third-order */
/* (i.e. cubic) interpolation. */

/* Also calculates the value of the independent variable T at which D is a */
/* minimum, with respect to the epoch of D1. */

/* N.B. The routine assumes that only one minimum is present in the interval H. */
/* === */
/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mce_min__(doublereal *d0, doublereal *d1, doublereal *
	d0t, doublereal *d1t, doublereal *h__, doublereal *d2min, doublereal *
	tmin)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal a, b, c__, tau, temp;



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    if (*d0t * *h__ > 0. || *d1t * *h__ < 0.) {
	if (*d0 <= *d1) {
	    *d2min = *d0;
	    *tmin = -(*h__);
	} else {
	    *d2min = *d1;
	    *tmin = 0.;
	}
    } else {
	temp = (*d0 - *d1) * 6.;
	a = temp + *h__ * 3. * (*d0t + *d1t);
	b = temp + *h__ * 2. * (*d0t + *d1t * 2.);
	c__ = *h__ * *d1t;

/* Computing MAX */
	d__2 = b * b - a * 4. * c__;
	d__1 = sqrt((max(d__2,0.)));
	temp = (b + d_sign(&d__1, &b)) * -.5;
	if (temp == 0.) {
	    tau = 0.;
	} else {
	    tau = c__ / temp;
	}

/* Make sure TAU falls in the interval -1 < TAU < 0 */
	tau = min(tau,0.);
	tau = max(tau,-1.);

/* Calculate TMIN and D2MIN */
	*tmin = tau * *h__;
	temp = tau + 1.;
	*d2min = tau * tau * ((tau * 2. + 3.) * *d0 + temp * *h__ * *d0t) + 
		temp * temp * ((1. - tau * 2.) * *d1 + tau * *h__ * *d1t);

/* Make sure D2MIN is not negative */
	*d2min = max(*d2min,0.);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mce_min__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*    MCE_SNIF.FOR    (ErikSoft   3 October 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Given initial and final coordinates and velocities X and V, and a timestep */
/* H, the routine estimates which objects were involved in a close encounter */
/* during the timestep. The routine examines all objects with indices I >= I0. */

/* Returns an array CE, which for each object is: */
/*                           0 if it will undergo no encounters */
/*                           2 if it will pass within RCRIT of a Big body */

/* Also returns arrays ICE and JCE, containing the indices of each pair of */
/* objects estimated to have undergone an encounter. */

/* N.B. All coordinates must be with respect to the central body!!!! */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mce_snif__(doublereal *h__, integer *i0, integer *nbod, 
	integer *nbig, doublereal *x0, doublereal *v0, doublereal *x1, 
	doublereal *v1, doublereal *rcrit, integer *ce, integer *nce, integer 
	*ice, integer *jce)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal d0, d1, rc, rc2, d0t, d1t, du0, dv0, dx0, dy0, dz0, dw0,
	     dx1, dy1, dz1, du1, dv1, dw1, temp, tmin, xmin[2000], ymin[2000],
	     xmax[2000], ymax[2000], d2min;
    extern /* Subroutine */ int mce_min__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), mce_box__(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);



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
    --ce;
    --rcrit;
    v1 -= 4;
    x1 -= 4;
    v0 -= 4;
    x0 -= 4;
    --ice;
    --jce;

    /* Function Body */
    if (*i0 <= 0) {
	*i0 = 2;
    }
    *nce = 0;
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	ce[j] = 0;
    }

/* Calculate maximum and minimum values of x and y coordinates */
    mce_box__(nbod, h__, &x0[4], &v0[4], &x1[4], &v1[4], xmin, xmax, ymin, 
	    ymax);

/* Adjust values for the Big bodies by symplectic close-encounter distance */
    i__1 = *nbig;
    for (j = *i0; j <= i__1; ++j) {
	xmin[j - 1] -= rcrit[j];
	xmax[j - 1] += rcrit[j];
	ymin[j - 1] -= rcrit[j];
	ymax[j - 1] += rcrit[j];
    }

/* Identify pairs whose X-Y boxes overlap, and calculate minimum separation */
    i__1 = *nbig;
    for (i__ = *i0; i__ <= i__1; ++i__) {
	i__2 = *nbod;
	for (j = i__ + 1; j <= i__2; ++j) {
	    if (xmax[i__ - 1] >= xmin[j - 1] && xmax[j - 1] >= xmin[i__ - 1] 
		    && ymax[i__ - 1] >= ymin[j - 1] && ymax[j - 1] >= ymin[
		    i__ - 1]) {

/* Determine the maximum separation that would qualify as an encounter */
/* Computing MAX */
		d__1 = rcrit[i__], d__2 = rcrit[j];
		rc = max(d__1,d__2);
		rc2 = rc * rc;

/* Calculate initial and final separations */
		dx0 = x0[i__ * 3 + 1] - x0[j * 3 + 1];
		dy0 = x0[i__ * 3 + 2] - x0[j * 3 + 2];
		dz0 = x0[i__ * 3 + 3] - x0[j * 3 + 3];
		dx1 = x1[i__ * 3 + 1] - x1[j * 3 + 1];
		dy1 = x1[i__ * 3 + 2] - x1[j * 3 + 2];
		dz1 = x1[i__ * 3 + 3] - x1[j * 3 + 3];
		d0 = dx0 * dx0 + dy0 * dy0 + dz0 * dz0;
		d1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;

/* Check for a possible minimum in between */
		du0 = v0[i__ * 3 + 1] - v0[j * 3 + 1];
		dv0 = v0[i__ * 3 + 2] - v0[j * 3 + 2];
		dw0 = v0[i__ * 3 + 3] - v0[j * 3 + 3];
		du1 = v1[i__ * 3 + 1] - v1[j * 3 + 1];
		dv1 = v1[i__ * 3 + 2] - v1[j * 3 + 2];
		dw1 = v1[i__ * 3 + 3] - v1[j * 3 + 3];
		d0t = (dx0 * du0 + dy0 * dv0 + dz0 * dw0) * 2.;
		d1t = (dx1 * du1 + dy1 * dv1 + dz1 * dw1) * 2.;

/* If separation derivative changes sign, find the minimum separation */
		d2min = 9.9e29;
		if (d0t * *h__ <= 0. && d1t * *h__ >= 0.) {
		    mce_min__(&d0, &d1, &d0t, &d1t, h__, &d2min, &tmin);
		}

/* If minimum separation is small enough, flag this as a possible encounter */
/* Computing MIN */
		d__1 = min(d0,d1);
		temp = min(d__1,d2min);
		if (temp <= rc2) {
		    ce[i__] = 2;
		    ce[j] = 2;
		    ++(*nce);
		    ice[*nce] = i__;
		    jce[*nce] = j;
		}
	    }
	}
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mce_snif__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*    MCE_STAT.FOR    (ErikSoft   1 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Returns details of all close-encounter minima involving at least one Big */
/* body during a timestep. The routine estimates minima using the initial */
/* and final coordinates X(0),X(1) and velocities V(0),V(1) of the step, and */
/* the stepsize H. */

/*  ICLO, JCLO contain the indices of the objects */
/*  DCLO is their minimum separation */
/*  TCLO is the time of closest approach relative to current time */

/* The routine also checks for collisions/near misses given the physical radii */
/* RPHYS, and returns the time THIT of the collision/near miss closest to the */
/* start of the timestep, and the identities IHIT and JHIT of the objects */
/* involved. */

/*  NHIT = +1 implies a collision */
/*         -1    "    a near miss */

/* N.B. All coordinates & velocities must be with respect to the central body!! */
/* === */
/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mce_stat__(doublereal *time, doublereal *h__, doublereal 
	*rcen, integer *nbod, integer *nbig, doublereal *m, doublereal *x0, 
	doublereal *v0, doublereal *x1, doublereal *v1, doublereal *rce, 
	doublereal *rphys, integer *nclo, integer *iclo, integer *jclo, 
	doublereal *dclo, doublereal *tclo, doublereal *ixvclo, doublereal *
	jxvclo, integer *nhit, integer *ihit, integer *jhit, integer *chit, 
	doublereal *dhit, doublereal *thit, doublereal *thit1, integer *
	nowflag, integer *stat, char *outfile, char *mem, integer *lmem, 
	ftnlen outfile_len, ftnlen mem_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void), f_clos(cllist *);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal d0, d1, hm1, d0t, d1t, du0, dv0, dx0, dy0, dz0, dw0, 
	    dx1, dy1, dz1, du1, dv1, dw1, d2ce, tmp0, tmp1, temp, tmin, xmin[
	    2000], xmax[2000], ymin[2000], ymax[2000], d2min, d2hit, d2near;
    extern /* Subroutine */ int mce_min__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), mce_box__(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___263 = { 0, 23, 0, "(/,2a,/,a)", 0 };




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
    --stat;
    --rphys;
    --rce;
    v1 -= 4;
    x1 -= 4;
    v0 -= 4;
    x0 -= 4;
    --m;
    --iclo;
    --jclo;
    --dclo;
    --tclo;
    ixvclo -= 7;
    jxvclo -= 7;
    --ihit;
    --jhit;
    --chit;
    --dhit;
    --thit;
    mem -= 80;
    --lmem;

    /* Function Body */
    *nhit = 0;
    *thit1 = d_sign(&c_b121, h__);
    hm1 = 1. / *h__;

/* Calculate maximum and minimum values of x and y coords for each object */
    mce_box__(nbod, h__, &x0[4], &v0[4], &x1[4], &v1[4], xmin, xmax, ymin, 
	    ymax);

/* Adjust values by the maximum close-encounter radius plus a fudge factor */
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	temp = rce[j] * 1.2;
	xmin[j - 1] -= temp;
	xmax[j - 1] += temp;
	ymin[j - 1] -= temp;
	ymax[j - 1] += temp;
    }

/* Check for close encounters between each pair of objects */
    i__1 = *nbig;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = *nbod;
	for (j = i__ + 1; j <= i__2; ++j) {
	    if (xmax[i__ - 1] >= xmin[j - 1] && xmax[j - 1] >= xmin[i__ - 1] 
		    && ymax[i__ - 1] >= ymin[j - 1] && ymax[j - 1] >= ymin[
		    i__ - 1] && stat[i__] >= 0 && stat[j] >= 0) {

/* If the X-Y boxes for this pair overlap, check circumstances more closely */
		dx0 = x0[i__ * 3 + 1] - x0[j * 3 + 1];
		dy0 = x0[i__ * 3 + 2] - x0[j * 3 + 2];
		dz0 = x0[i__ * 3 + 3] - x0[j * 3 + 3];
		du0 = v0[i__ * 3 + 1] - v0[j * 3 + 1];
		dv0 = v0[i__ * 3 + 2] - v0[j * 3 + 2];
		dw0 = v0[i__ * 3 + 3] - v0[j * 3 + 3];
		d0t = (dx0 * du0 + dy0 * dv0 + dz0 * dw0) * 2.;

		dx1 = x1[i__ * 3 + 1] - x1[j * 3 + 1];
		dy1 = x1[i__ * 3 + 2] - x1[j * 3 + 2];
		dz1 = x1[i__ * 3 + 3] - x1[j * 3 + 3];
		du1 = v1[i__ * 3 + 1] - v1[j * 3 + 1];
		dv1 = v1[i__ * 3 + 2] - v1[j * 3 + 2];
		dw1 = v1[i__ * 3 + 3] - v1[j * 3 + 3];
		d1t = (dx1 * du1 + dy1 * dv1 + dz1 * dw1) * 2.;

/* Estimate minimum separation during the time interval, using interpolation */
		d0 = dx0 * dx0 + dy0 * dy0 + dz0 * dz0;
		d1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
		mce_min__(&d0, &d1, &d0t, &d1t, h__, &d2min, &tmin);
/* Computing MAX */
		d__1 = rce[i__], d__2 = rce[j];
		d2ce = max(d__1,d__2);
		d2hit = rphys[i__] + rphys[j];
		d2ce *= d2ce;
		d2hit *= d2hit;
		d2near = d2hit * 4.;

/* If the minimum separation qualifies as an encounter or if a collision */
/* is in progress, store details */
		if (d2min <= d2ce && d0t * *h__ <= 0. && d1t * *h__ >= 0. || 
			d2min <= d2hit) {
		    ++(*nclo);
		    if (*nclo > 50) {
L230:
			o__1.oerr = 1;
			o__1.ounit = 23;
			o__1.ofnmlen = 80;
			o__1.ofnm = outfile;
			o__1.orl = 0;
			o__1.osta = "old";
			o__1.oacc = "append";
			o__1.ofm = 0;
			o__1.oblnk = 0;
			i__3 = f_open(&o__1);
			if (i__3 != 0) {
			    goto L230;
			}
			s_wsfe(&io___263);
			do_fio(&c__1, mem + 9680, lmem[121]);
			do_fio(&c__1, mem + 10560, lmem[132]);
			do_fio(&c__1, mem + 6560, lmem[82]);
			e_wsfe();
			cl__1.cerr = 0;
			cl__1.cunit = 23;
			cl__1.csta = 0;
			f_clos(&cl__1);
		    } else {
			tclo[*nclo] = tmin + *time;
			dclo[*nclo] = sqrt((max(0.,d2min)));
			iclo[*nclo] = i__;
			jclo[*nclo] = j;

/* Make sure the more massive body is listed first */
			if (m[j] > m[i__] && j <= *nbig) {
			    iclo[*nclo] = j;
			    jclo[*nclo] = i__;
			}

/* Make linear interpolation to get coordinates at time of closest approach */
			tmp0 = tmin * hm1 + 1.;
			tmp1 = -tmin * hm1;
			ixvclo[*nclo * 6 + 1] = tmp0 * x0[i__ * 3 + 1] + tmp1 
				* x1[i__ * 3 + 1];
			ixvclo[*nclo * 6 + 2] = tmp0 * x0[i__ * 3 + 2] + tmp1 
				* x1[i__ * 3 + 2];
			ixvclo[*nclo * 6 + 3] = tmp0 * x0[i__ * 3 + 3] + tmp1 
				* x1[i__ * 3 + 3];
			ixvclo[*nclo * 6 + 4] = tmp0 * v0[i__ * 3 + 1] + tmp1 
				* v1[i__ * 3 + 1];
			ixvclo[*nclo * 6 + 5] = tmp0 * v0[i__ * 3 + 2] + tmp1 
				* v1[i__ * 3 + 2];
			ixvclo[*nclo * 6 + 6] = tmp0 * v0[i__ * 3 + 3] + tmp1 
				* v1[i__ * 3 + 3];
			jxvclo[*nclo * 6 + 1] = tmp0 * x0[j * 3 + 1] + tmp1 * 
				x1[j * 3 + 1];
			jxvclo[*nclo * 6 + 2] = tmp0 * x0[j * 3 + 2] + tmp1 * 
				x1[j * 3 + 2];
			jxvclo[*nclo * 6 + 3] = tmp0 * x0[j * 3 + 3] + tmp1 * 
				x1[j * 3 + 3];
			jxvclo[*nclo * 6 + 4] = tmp0 * v0[j * 3 + 1] + tmp1 * 
				v1[j * 3 + 1];
			jxvclo[*nclo * 6 + 5] = tmp0 * v0[j * 3 + 2] + tmp1 * 
				v1[j * 3 + 2];
			jxvclo[*nclo * 6 + 6] = tmp0 * v0[j * 3 + 3] + tmp1 * 
				v1[j * 3 + 3];
		    }
		}

/* Check for a near miss or collision */
		if (d2min <= d2near) {
		    ++(*nhit);
		    ihit[*nhit] = i__;
		    jhit[*nhit] = j;
		    thit[*nhit] = tmin + *time;
		    dhit[*nhit] = sqrt(d2min);
		    chit[*nhit] = -1;
		    if (d2min <= d2hit) {
			chit[*nhit] = 1;
		    }

/* Make sure the more massive body is listed first */
		    if (m[jhit[*nhit]] > m[ihit[*nhit]] && j <= *nbig) {
			ihit[*nhit] = j;
			jhit[*nhit] = i__;
		    }

/* Is this the collision closest to the start of the time step? */
		    if ((tmin - *thit1) * *h__ < 0.) {
			*thit1 = tmin;
			*nowflag = 0;
			if (d1 <= d2hit) {
			    *nowflag = 1;
			}
		    }
		}
	    }

/* Move on to the next pair of objects */
	}
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mce_stat__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_ACSH.FOR    (ErikSoft  2 March 1999) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates inverse hyperbolic cosine of an angle X (in radians). */

/* ------------------------------------------------------------------------------ */

doublereal mco_acsh__(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal);



/* Input/Output */

/* ------------------------------------------------------------------------------ */

    if (*x >= 1.) {
	ret_val = log(*x + sqrt(*x * *x - 1.));
    } else {
	ret_val = 0.;
    }

/* ------------------------------------------------------------------------------ */

    return ret_val;
} /* mco_acsh__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_B2H.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Converts barycentric coordinates to coordinates with respect to the central */
/* body. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_b2h__(doublereal *time, doublereal *jcen, integer *
	nbod, integer *nbig, doublereal *h__, doublereal *m, doublereal *x, 
	doublereal *v, doublereal *xh, doublereal *vh, doublereal *ngf, 
	integer *ngflag, integer *opt)
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
    ngf -= 5;
    vh -= 4;
    xh -= 4;
    v -= 4;
    x -= 4;
    --m;
    --opt;

    /* Function Body */
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	xh[j * 3 + 1] = x[j * 3 + 1] - x[4];
	xh[j * 3 + 2] = x[j * 3 + 2] - x[5];
	xh[j * 3 + 3] = x[j * 3 + 3] - x[6];
	vh[j * 3 + 1] = v[j * 3 + 1] - v[4];
	vh[j * 3 + 2] = v[j * 3 + 2] - v[5];
	vh[j * 3 + 3] = v[j * 3 + 3] - v[6];
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_b2h__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_DH2H.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Converts democratic heliocentric coordinates to coordinates with respect */
/* to the central body. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_dh2h__(doublereal *time, doublereal *jcen, integer *
	nbod, integer *nbig, doublereal *h__, doublereal *m, doublereal *x, 
	doublereal *v, doublereal *xh, doublereal *vh, doublereal *ngf, 
	integer *ngflag, integer *opt)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal temp, mvsum[3];



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    --jcen;
    ngf -= 5;
    vh -= 4;
    xh -= 4;
    v -= 4;
    x -= 4;
    --m;
    --opt;

    /* Function Body */
    mvsum[0] = 0.;
    mvsum[1] = 0.;
    mvsum[2] = 0.;

    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	xh[j * 3 + 1] = x[j * 3 + 1];
	xh[j * 3 + 2] = x[j * 3 + 2];
	xh[j * 3 + 3] = x[j * 3 + 3];
	mvsum[0] += m[j] * v[j * 3 + 1];
	mvsum[1] += m[j] * v[j * 3 + 2];
	mvsum[2] += m[j] * v[j * 3 + 3];
    }

    temp = 1. / m[1];
    mvsum[0] = temp * mvsum[0];
    mvsum[1] = temp * mvsum[1];
    mvsum[2] = temp * mvsum[2];

    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	vh[j * 3 + 1] = v[j * 3 + 1] + mvsum[0];
	vh[j * 3 + 2] = v[j * 3 + 2] + mvsum[1];
	vh[j * 3 + 3] = v[j * 3 + 3] + mvsum[2];
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_dh2h__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_IDEN.FOR    (ErikSoft   2 November 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Makes a new copy of a set of coordinates. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_iden__(doublereal *time, doublereal *jcen, integer *
	nbod, integer *nbig, doublereal *h__, doublereal *m, doublereal *xh, 
	doublereal *vh, doublereal *x, doublereal *v, doublereal *ngf, 
	integer *ngflag, integer *opt)
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
    ngf -= 5;
    v -= 4;
    x -= 4;
    vh -= 4;
    xh -= 4;
    --m;
    --opt;

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

/*      MCO_MVS2H.FOR    (ErikSoft   28 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Applies a symplectic corrector, which converts coordinates for a second- */
/* order mixed-variable symplectic integrator to ones with respect to the */
/* central body. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_mvs2h__(doublereal *time, doublereal *jcen, integer *
	nbod, integer *nbig, doublereal *h__, doublereal *m, doublereal *x, 
	doublereal *v, doublereal *xh, doublereal *vh, doublereal *ngf, 
	integer *ngflag, integer *opt)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int mco_iden__(doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), mfo_user__(doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    ;
    static doublereal a[6000]	/* was [3][2000] */;
    static integer j, k;
    extern /* Subroutine */ int drift_one__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal ha[2], hb[2], gm[2000], vj[6000]	/* was [3][2000] */, 
	    xj[6000]	/* was [3][2000] */, rt10, angf[6000]	/* was [3][
	    2000] */, ausr[6000]	/* was [3][2000] */;
    static integer stat[2000], iflag;
    static doublereal msofar;
    extern /* Subroutine */ int mco_h2j__(doublereal *, doublereal *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), mco_j2h__(doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    mfo_ngf__(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal minside;
    extern /* Subroutine */ int mfo_mvs__(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;



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
    --jcen;
    ngf -= 5;
    vh -= 4;
    xh -= 4;
    v -= 4;
    x -= 4;
    --m;
    --opt;

    /* Function Body */
    rt10 = sqrt(10.);
    ha[0] = *h__ * rt10 * 3. / 10.;
    hb[0] = -(*h__) * rt10 / 72.;
    ha[1] = *h__ * rt10 / 5.;
    hb[1] = *h__ * rt10 / 24.;
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	angf[j * 3 - 3] = 0.;
	angf[j * 3 - 2] = 0.;
	angf[j * 3 - 1] = 0.;
	ausr[j * 3 - 3] = 0.;
	ausr[j * 3 - 2] = 0.;
	ausr[j * 3 - 1] = 0.;
    }
    mco_iden__(time, &jcen[1], nbod, nbig, h__, &m[1], &x[4], &v[4], &xh[4], &
	    vh[4], &ngf[5], ngflag, &opt[1]);

/* Calculate effective central masses for Kepler drifts */
    minside = m[1];
    i__1 = *nbig;
    for (j = 2; j <= i__1; ++j) {
	msofar = minside + m[j];
	gm[j - 1] = m[1] * msofar / minside;
	minside = msofar;
    }

/* Two step corrector */
    for (k = 1; k <= 2; ++k) {

/* Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies) */
	mco_h2j__(time, &jcen[1], nbig, nbig, h__, &m[1], &xh[4], &vh[4], xj, 
		vj, &ngf[5], ngflag, &opt[1]);
	i__1 = *nbig;
	for (j = 2; j <= i__1; ++j) {
	    iflag = 0;
	    drift_one__(&gm[j - 1], &xj[j * 3 - 3], &xj[j * 3 - 2], &xj[j * 3 
		    - 1], &vj[j * 3 - 3], &vj[j * 3 - 2], &vj[j * 3 - 1], &ha[
		    k - 1], &iflag);
	}
	i__1 = *nbod;
	for (j = *nbig + 1; j <= i__1; ++j) {
	    iflag = 0;
	    drift_one__(&m[1], &xh[j * 3 + 1], &xh[j * 3 + 2], &xh[j * 3 + 3],
		     &vh[j * 3 + 1], &vh[j * 3 + 2], &vh[j * 3 + 3], &ha[k - 
		    1], &iflag);
	}

/* Advance Interaction Hamiltonian */
	mco_j2h__(time, &jcen[1], nbig, nbig, h__, &m[1], xj, vj, &xh[4], &vh[
		4], &ngf[5], ngflag, &opt[1]);
	mfo_mvs__(&jcen[1], nbod, nbig, &m[1], &xh[4], xj, a, stat);

/* If required, apply non-gravitational and user-defined forces */
	if (opt[8] == 1) {
	    mfo_user__(time, &jcen[1], nbod, nbig, &m[1], &xh[4], &vh[4], 
		    ausr);
	}
	if (*ngflag == 1 || *ngflag == 3) {
	    mfo_ngf__(nbod, &xh[4], &vh[4], angf, &ngf[5]);
	}

	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    vh[j * 3 + 1] -= hb[k - 1] * (angf[j * 3 - 3] + ausr[j * 3 - 3] + 
		    a[j * 3 - 3]);
	    vh[j * 3 + 2] -= hb[k - 1] * (angf[j * 3 - 2] + ausr[j * 3 - 2] + 
		    a[j * 3 - 2]);
	    vh[j * 3 + 3] -= hb[k - 1] * (angf[j * 3 - 1] + ausr[j * 3 - 1] + 
		    a[j * 3 - 1]);
	}

/* Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies) */
	mco_h2j__(time, &jcen[1], nbig, nbig, h__, &m[1], &xh[4], &vh[4], xj, 
		vj, &ngf[5], ngflag, &opt[1]);
	i__1 = *nbig;
	for (j = 2; j <= i__1; ++j) {
	    iflag = 0;
	    d__1 = ha[k - 1] * -2.;
	    drift_one__(&gm[j - 1], &xj[j * 3 - 3], &xj[j * 3 - 2], &xj[j * 3 
		    - 1], &vj[j * 3 - 3], &vj[j * 3 - 2], &vj[j * 3 - 1], &
		    d__1, &iflag);
	}
	i__1 = *nbod;
	for (j = *nbig + 1; j <= i__1; ++j) {
	    iflag = 0;
	    d__1 = ha[k - 1] * -2.;
	    drift_one__(&m[1], &xh[j * 3 + 1], &xh[j * 3 + 2], &xh[j * 3 + 3],
		     &vh[j * 3 + 1], &vh[j * 3 + 2], &vh[j * 3 + 3], &d__1, &
		    iflag);
	}

/* Advance Interaction Hamiltonian */
	mco_j2h__(time, &jcen[1], nbig, nbig, h__, &m[1], xj, vj, &xh[4], &vh[
		4], &ngf[5], ngflag, &opt[1]);
	mfo_mvs__(&jcen[1], nbod, nbig, &m[1], &xh[4], xj, a, stat);

/* If required, apply non-gravitational and user-defined forces */
	if (opt[8] == 1) {
	    mfo_user__(time, &jcen[1], nbod, nbig, &m[1], &xh[4], &vh[4], 
		    ausr);
	}
	if (*ngflag == 1 || *ngflag == 3) {
	    mfo_ngf__(nbod, &xh[4], &vh[4], angf, &ngf[5]);
	}

	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    vh[j * 3 + 1] += hb[k - 1] * (angf[j * 3 - 3] + ausr[j * 3 - 3] + 
		    a[j * 3 - 3]);
	    vh[j * 3 + 2] += hb[k - 1] * (angf[j * 3 - 2] + ausr[j * 3 - 2] + 
		    a[j * 3 - 2]);
	    vh[j * 3 + 3] += hb[k - 1] * (angf[j * 3 - 1] + ausr[j * 3 - 1] + 
		    a[j * 3 - 1]);
	}

/* Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies) */
	mco_h2j__(time, &jcen[1], nbig, nbig, h__, &m[1], &xh[4], &vh[4], xj, 
		vj, &ngf[5], ngflag, &opt[1]);
	i__1 = *nbig;
	for (j = 2; j <= i__1; ++j) {
	    iflag = 0;
	    drift_one__(&gm[j - 1], &xj[j * 3 - 3], &xj[j * 3 - 2], &xj[j * 3 
		    - 1], &vj[j * 3 - 3], &vj[j * 3 - 2], &vj[j * 3 - 1], &ha[
		    k - 1], &iflag);
	}
	i__1 = *nbod;
	for (j = *nbig + 1; j <= i__1; ++j) {
	    iflag = 0;
	    drift_one__(&m[1], &xh[j * 3 + 1], &xh[j * 3 + 2], &xh[j * 3 + 3],
		     &vh[j * 3 + 1], &vh[j * 3 + 2], &vh[j * 3 + 3], &ha[k - 
		    1], &iflag);
	}
	mco_j2h__(time, &jcen[1], nbig, nbig, h__, &m[1], xj, vj, &xh[4], &vh[
		4], &ngf[5], ngflag, &opt[1]);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_mvs2h__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_EL2X.FOR    (ErikSoft  7 July 1999) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates Cartesian coordinates and velocities given Keplerian orbital */
/* elements (for elliptical, parabolic or hyperbolic orbits). */

/* Based on a routine from Levison and Duncan's SWIFT integrator. */

/*  gm = grav const * (central + secondary mass) */
/*  q = perihelion distance */
/*  e = eccentricity */
/*  i = inclination                 ) */
/*  p = longitude of perihelion !!! )   in */
/*  n = longitude of ascending node ) radians */
/*  l = mean anomaly                ) */

/*  x,y,z = Cartesian positions  ( units the same as a ) */
/*  u,v,w =     "     velocities ( units the same as sqrt(gm/a) ) */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_el2x__(doublereal *gm, doublereal *q, doublereal *e, 
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
	temp = sqrt(*gm / a) / (1. - *e * ce);
	z3 = -se * temp;
	z4 = romes * ce * temp;
    } else {
/* Parabola */
	if (*e == 1.) {
	    ce = orbel_zget__(l);
	    z1 = *q * (1. - ce * ce);
	    z2 = *q * 2. * ce;
	    z4 = sqrt(*gm * 2. / *q) / (ce * ce + 1.);
	    z3 = -ce * z4;
	} else {
/* Hyperbola */
	    romes = sqrt(*e * *e - 1.);
	    temp = orbel_fhybrid__(e, l);
	    mco_sinh__(&temp, &se, &ce);
	    z1 = a * (ce - *e);
	    z2 = -a * romes * se;
	    temp = sqrt(*gm / abs(a)) / (*e * ce - 1.);
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

/*      MCO_H2B.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Converts coordinates with respect to the central body to barycentric */
/* coordinates. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_h2b__(doublereal *time, doublereal *jcen, integer *
	nbod, integer *nbig, doublereal *h__, doublereal *m, doublereal *xh, 
	doublereal *vh, doublereal *x, doublereal *v, doublereal *ngf, 
	integer *ngflag, integer *opt)
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
    ngf -= 5;
    v -= 4;
    x -= 4;
    vh -= 4;
    xh -= 4;
    --m;
    --opt;

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

/*      MCO_H2MVS.FOR    (ErikSoft   28 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Applies an inverse symplectic corrector, which converts coordinates with */
/* respect to the central body to integrator coordinates for a second-order */
/* mixed-variable symplectic integrator. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_h2mvs__(doublereal *time, doublereal *jcen, integer *
	nbod, integer *nbig, doublereal *h__, doublereal *m, doublereal *xh, 
	doublereal *vh, doublereal *x, doublereal *v, doublereal *ngf, 
	integer *ngflag, integer *opt)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int mco_iden__(doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), mfo_user__(doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    ;
    static doublereal a[6000]	/* was [3][2000] */;
    static integer j, k;
    extern /* Subroutine */ int drift_one__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal ha[2], hb[2], gm[2000], vj[6000]	/* was [3][2000] */, 
	    xj[6000]	/* was [3][2000] */, rt10, angf[6000]	/* was [3][
	    2000] */, ausr[6000]	/* was [3][2000] */;
    static integer stat[2000], iflag;
    static doublereal msofar;
    extern /* Subroutine */ int mco_h2j__(doublereal *, doublereal *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), mco_j2h__(doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    mfo_ngf__(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal minside;
    extern /* Subroutine */ int mfo_mvs__(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;



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
    --jcen;
    ngf -= 5;
    v -= 4;
    x -= 4;
    vh -= 4;
    xh -= 4;
    --m;
    --opt;

    /* Function Body */
    rt10 = sqrt(10.);
    ha[0] = -(*h__) * rt10 / 5.;
    hb[0] = -(*h__) * rt10 / 24.;
    ha[1] = -(*h__) * rt10 * 3. / 10.;
    hb[1] = *h__ * rt10 / 72.;
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	angf[j * 3 - 3] = 0.;
	angf[j * 3 - 2] = 0.;
	angf[j * 3 - 1] = 0.;
	ausr[j * 3 - 3] = 0.;
	ausr[j * 3 - 2] = 0.;
	ausr[j * 3 - 1] = 0.;
    }
    mco_iden__(time, &jcen[1], nbod, nbig, h__, &m[1], &xh[4], &vh[4], &x[4], 
	    &v[4], &ngf[5], ngflag, &opt[1]);

/* Calculate effective central masses for Kepler drifts */
    minside = m[1];
    i__1 = *nbig;
    for (j = 2; j <= i__1; ++j) {
	msofar = minside + m[j];
	gm[j - 1] = m[1] * msofar / minside;
	minside = msofar;
    }

    for (k = 1; k <= 2; ++k) {

/* Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies) */
	mco_h2j__(time, &jcen[1], nbig, nbig, h__, &m[1], &x[4], &v[4], xj, 
		vj, &ngf[5], ngflag, &opt[1]);
	i__1 = *nbig;
	for (j = 2; j <= i__1; ++j) {
	    iflag = 0;
	    drift_one__(&gm[j - 1], &xj[j * 3 - 3], &xj[j * 3 - 2], &xj[j * 3 
		    - 1], &vj[j * 3 - 3], &vj[j * 3 - 2], &vj[j * 3 - 1], &ha[
		    k - 1], &iflag);
	}
	i__1 = *nbod;
	for (j = *nbig + 1; j <= i__1; ++j) {
	    iflag = 0;
	    drift_one__(&m[1], &x[j * 3 + 1], &x[j * 3 + 2], &x[j * 3 + 3], &
		    v[j * 3 + 1], &v[j * 3 + 2], &v[j * 3 + 3], &ha[k - 1], &
		    iflag);
	}

/* Advance Interaction Hamiltonian */
	mco_j2h__(time, &jcen[1], nbig, nbig, h__, &m[1], xj, vj, &x[4], &v[4]
		, &ngf[5], ngflag, &opt[1]);
	mfo_mvs__(&jcen[1], nbod, nbig, &m[1], &x[4], xj, a, stat);

/* If required, apply non-gravitational and user-defined forces */
	if (opt[8] == 1) {
	    mfo_user__(time, &jcen[1], nbod, nbig, &m[1], &x[4], &v[4], ausr);
	}
	if (*ngflag == 1 || *ngflag == 3) {
	    mfo_ngf__(nbod, &x[4], &v[4], angf, &ngf[5]);
	}

	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    v[j * 3 + 1] += hb[k - 1] * (angf[j * 3 - 3] + ausr[j * 3 - 3] + 
		    a[j * 3 - 3]);
	    v[j * 3 + 2] += hb[k - 1] * (angf[j * 3 - 2] + ausr[j * 3 - 2] + 
		    a[j * 3 - 2]);
	    v[j * 3 + 3] += hb[k - 1] * (angf[j * 3 - 1] + ausr[j * 3 - 1] + 
		    a[j * 3 - 1]);
	}

/* Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies) */
	mco_h2j__(time, &jcen[1], nbig, nbig, h__, &m[1], &x[4], &v[4], xj, 
		vj, &ngf[5], ngflag, &opt[1]);
	i__1 = *nbig;
	for (j = 2; j <= i__1; ++j) {
	    iflag = 0;
	    d__1 = ha[k - 1] * -2.;
	    drift_one__(&gm[j - 1], &xj[j * 3 - 3], &xj[j * 3 - 2], &xj[j * 3 
		    - 1], &vj[j * 3 - 3], &vj[j * 3 - 2], &vj[j * 3 - 1], &
		    d__1, &iflag);
	}
	i__1 = *nbod;
	for (j = *nbig + 1; j <= i__1; ++j) {
	    iflag = 0;
	    d__1 = ha[k - 1] * -2.;
	    drift_one__(&m[1], &x[j * 3 + 1], &x[j * 3 + 2], &x[j * 3 + 3], &
		    v[j * 3 + 1], &v[j * 3 + 2], &v[j * 3 + 3], &d__1, &iflag)
		    ;
	}

/* Advance Interaction Hamiltonian */
	mco_j2h__(time, &jcen[1], nbig, nbig, h__, &m[1], xj, vj, &x[4], &v[4]
		, &ngf[5], ngflag, &opt[1]);
	mfo_mvs__(&jcen[1], nbod, nbig, &m[1], &x[4], xj, a, stat);

/* If required, apply non-gravitational and user-defined forces */
	if (opt[8] == 1) {
	    mfo_user__(time, &jcen[1], nbod, nbig, &m[1], &x[4], &v[4], ausr);
	}
	if (*ngflag == 1 || *ngflag == 3) {
	    mfo_ngf__(nbod, &x[4], &v[4], angf, &ngf[5]);
	}

	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    v[j * 3 + 1] -= hb[k - 1] * (angf[j * 3 - 3] + ausr[j * 3 - 3] + 
		    a[j * 3 - 3]);
	    v[j * 3 + 2] -= hb[k - 1] * (angf[j * 3 - 2] + ausr[j * 3 - 2] + 
		    a[j * 3 - 2]);
	    v[j * 3 + 3] -= hb[k - 1] * (angf[j * 3 - 1] + ausr[j * 3 - 1] + 
		    a[j * 3 - 1]);
	}

/* Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies) */
	mco_h2j__(time, &jcen[1], nbig, nbig, h__, &m[1], &x[4], &v[4], xj, 
		vj, &ngf[5], ngflag, &opt[1]);
	i__1 = *nbig;
	for (j = 2; j <= i__1; ++j) {
	    iflag = 0;
	    drift_one__(&gm[j - 1], &xj[j * 3 - 3], &xj[j * 3 - 2], &xj[j * 3 
		    - 1], &vj[j * 3 - 3], &vj[j * 3 - 2], &vj[j * 3 - 1], &ha[
		    k - 1], &iflag);
	}
	i__1 = *nbod;
	for (j = *nbig + 1; j <= i__1; ++j) {
	    iflag = 0;
	    drift_one__(&m[1], &x[j * 3 + 1], &x[j * 3 + 2], &x[j * 3 + 3], &
		    v[j * 3 + 1], &v[j * 3 + 2], &v[j * 3 + 3], &ha[k - 1], &
		    iflag);
	}
	mco_j2h__(time, &jcen[1], nbig, nbig, h__, &m[1], xj, vj, &x[4], &v[4]
		, &ngf[5], ngflag, &opt[1]);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_h2mvs__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_H2DH.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Convert coordinates with respect to the central body to democratic */
/* heliocentric coordinates. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_h2dh__(doublereal *time, doublereal *jcen, integer *
	nbod, integer *nbig, doublereal *h__, doublereal *m, doublereal *xh, 
	doublereal *vh, doublereal *x, doublereal *v, doublereal *ngf, 
	integer *ngflag, integer *opt)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal temp, mtot, mvsum[3];



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    --jcen;
    ngf -= 5;
    v -= 4;
    x -= 4;
    vh -= 4;
    xh -= 4;
    --m;
    --opt;

    /* Function Body */
    mtot = 0.;
    mvsum[0] = 0.;
    mvsum[1] = 0.;
    mvsum[2] = 0.;

    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	x[j * 3 + 1] = xh[j * 3 + 1];
	x[j * 3 + 2] = xh[j * 3 + 2];
	x[j * 3 + 3] = xh[j * 3 + 3];
	mtot += m[j];
	mvsum[0] += m[j] * vh[j * 3 + 1];
	mvsum[1] += m[j] * vh[j * 3 + 2];
	mvsum[2] += m[j] * vh[j * 3 + 3];
    }

    temp = 1. / (m[1] + mtot);

    mvsum[0] = temp * mvsum[0];
    mvsum[1] = temp * mvsum[1];
    mvsum[2] = temp * mvsum[2];

    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	v[j * 3 + 1] = vh[j * 3 + 1] - mvsum[0];
	v[j * 3 + 2] = vh[j * 3 + 2] - mvsum[1];
	v[j * 3 + 3] = vh[j * 3 + 3] - mvsum[2];
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_h2dh__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_H2J.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Converts coordinates with respect to the central body to Jacobi coordinates. */

/* N.B. The coordinates respect to the central body for the small bodies */
/* ===  are assumed to be equal to their Jacobi coordinates. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_h2j__(doublereal *time, doublereal *jcen, integer *
	nbod, integer *nbig, doublereal *h__, doublereal *m, doublereal *xh, 
	doublereal *vh, doublereal *x, doublereal *v, doublereal *ngf, 
	integer *ngflag, integer *opt)
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
    ngf -= 5;
    v -= 4;
    x -= 4;
    vh -= 4;
    xh -= 4;
    --m;
    --opt;

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

/*      MCO_J2H.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Converts Jacobi coordinates to coordinates with respect to the central */
/* body. */

/* N.B. The Jacobi coordinates of the small bodies are assumed to be equal */
/* ===  to their coordinates with respect to the central body. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_j2h__(doublereal *time, doublereal *jcen, integer *
	nbod, integer *nbig, doublereal *h__, doublereal *m, doublereal *x, 
	doublereal *v, doublereal *xh, doublereal *vh, doublereal *ngf, 
	integer *ngflag, integer *opt)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal mu, mv, mw, mx, my, mz, temp, mtot;



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    --jcen;
    ngf -= 5;
    vh -= 4;
    xh -= 4;
    v -= 4;
    x -= 4;
    --m;
    --opt;

    /* Function Body */
    xh[7] = x[7];
    xh[8] = x[8];
    xh[9] = x[9];
    vh[7] = v[7];
    vh[8] = v[8];
    vh[9] = v[9];
    mtot = m[2];
    temp = m[2] / (mtot + m[1]);
    mx = temp * x[7];
    my = temp * x[8];
    mz = temp * x[9];
    mu = temp * v[7];
    mv = temp * v[8];
    mw = temp * v[9];

    i__1 = *nbig - 1;
    for (j = 3; j <= i__1; ++j) {
	xh[j * 3 + 1] = x[j * 3 + 1] + mx;
	xh[j * 3 + 2] = x[j * 3 + 2] + my;
	xh[j * 3 + 3] = x[j * 3 + 3] + mz;
	vh[j * 3 + 1] = v[j * 3 + 1] + mu;
	vh[j * 3 + 2] = v[j * 3 + 2] + mv;
	vh[j * 3 + 3] = v[j * 3 + 3] + mw;
	mtot += m[j];
	temp = m[j] / (mtot + m[1]);
	mx += temp * x[j * 3 + 1];
	my += temp * x[j * 3 + 2];
	mz += temp * x[j * 3 + 3];
	mu += temp * v[j * 3 + 1];
	mv += temp * v[j * 3 + 2];
	mw += temp * v[j * 3 + 3];
    }

    if (*nbig > 2) {
	xh[*nbig * 3 + 1] = x[*nbig * 3 + 1] + mx;
	xh[*nbig * 3 + 2] = x[*nbig * 3 + 2] + my;
	xh[*nbig * 3 + 3] = x[*nbig * 3 + 3] + mz;
	vh[*nbig * 3 + 1] = v[*nbig * 3 + 1] + mu;
	vh[*nbig * 3 + 2] = v[*nbig * 3 + 2] + mv;
	vh[*nbig * 3 + 3] = v[*nbig * 3 + 3] + mw;
    }

    i__1 = *nbod;
    for (j = *nbig + 1; j <= i__1; ++j) {
	xh[j * 3 + 1] = x[j * 3 + 1];
	xh[j * 3 + 2] = x[j * 3 + 2];
	xh[j * 3 + 3] = x[j * 3 + 3];
	vh[j * 3 + 1] = v[j * 3 + 1];
	vh[j * 3 + 2] = v[j * 3 + 2];
	vh[j * 3 + 3] = v[j * 3 + 3];
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_j2h__ */


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

/*      MCO_X2A.FOR    (ErikSoft   4 October 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates an object's orbital semi-major axis given its Cartesian coords. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_x2a__(doublereal *gm, doublereal *x, doublereal *y, 
	doublereal *z__, doublereal *u, doublereal *v, doublereal *w, 
	doublereal *a, doublereal *r__, doublereal *v2)
{
    /* Builtin functions */
    double sqrt(doublereal);



/* Input/Output */

/* ------------------------------------------------------------------------------ */

    *r__ = sqrt(*x * *x + *y * *y + *z__ * *z__);
    *v2 = *u * *u + *v * *v + *w * *w;
    *a = *gm * *r__ / (*gm * 2. - *r__ * *v2);

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_x2a__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_X2OV.FOR    (ErikSoft   20 February 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates output variables for an object given its coordinates and */
/* velocities. The output variables are: */
/*  r = the radial distance */
/*  theta = polar angle */
/*  phi = azimuthal angle */
/*  fv = 1 / [1 + 2(ke/be)^2], where be and ke are the object's binding and */
/*                             kinetic energies. (Note that 0 < fv < 1). */
/*  vtheta = polar angle of velocity vector */
/*  vphi = azimuthal angle of the velocity vector */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mco_x2ov__(doublereal *rcen, doublereal *rmax, 
	doublereal *mcen, doublereal *m, doublereal *x, doublereal *y, 
	doublereal *z__, doublereal *u, doublereal *v, doublereal *w, 
	doublereal *fr, doublereal *theta, doublereal *phi, doublereal *fv, 
	doublereal *vtheta, doublereal *vphi)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_lg10(doublereal *), acos(doublereal), d_mod(
	    doublereal *, doublereal *), atan2(doublereal, doublereal);

    /* Local variables */
    static doublereal r__, v1, v2, be, ke, temp;



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

    r__ = sqrt(*x * *x + *y * *y + *z__ * *z__);
    v2 = *u * *u + *v * *v + *w * *w;
    v1 = sqrt(v2);
    be = (*mcen + *m) / r__;
    ke = v2 * .5;

/* Computing MIN */
    d__2 = max(r__,*rcen);
    d__1 = min(d__2,*rmax) / *rcen;
    *fr = d_lg10(&d__1);
    temp = ke / be;
    *fv = 1. / (temp * 2. * temp + 1.);

    d__1 = acos(*z__ / r__) + 6.2831853071795862;
    *theta = d_mod(&d__1, &c_b58);
    d__1 = acos(*w / v1) + 6.2831853071795862;
    *vtheta = d_mod(&d__1, &c_b58);
    d__1 = atan2(*y, *x) + 6.2831853071795862;
    *phi = d_mod(&d__1, &c_b58);
    d__1 = atan2(*v, *u) + 6.2831853071795862;
    *vphi = d_mod(&d__1, &c_b58);

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_x2ov__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MCO_X2EL.FOR    (ErikSoft  23 January 2001) */

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
		ce = d_sign(&c_b150, &ce);
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
	    cf = d_sign(&c_b150, &cf);
	}
	f = acos(cf);
	if (rv < 0.) {
	    f = 6.2831853071795862 - f;
	}
	*p = true__ - f;
	d__1 = *p + 6.2831853071795862 + 6.2831853071795862;
	*p = d_mod(&d__1, &c_b58);
    }

    if (*l < 0.) {
	*l += 6.2831853071795862;
    }
    if (*l > 6.2831853071795862) {
	*l = d_mod(l, &c_b58);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mco_x2el__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MDT_BS1.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Integrates NBOD bodies (of which NBIG are Big) for one timestep H0 */
/* using the Bulirsch-Stoer method. The accelerations are calculated using the */
/* subroutine FORCE. The accuracy of the step is approximately determined */
/* by the tolerance parameter TOL. */

/* N.B. Input/output must be in coordinates with respect to the central body. */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mdt_bs1__(doublereal *time, doublereal *h0, doublereal *
	hdid, doublereal *tol, doublereal *jcen, integer *nbod, integer *nbig,
	 doublereal *mass, doublereal *x0, doublereal *v0, doublereal *s, 
	doublereal *rphys, doublereal *rcrit, doublereal *ngf, integer *stat, 
	integer *dtflag, integer *ngflag, integer *opt, integer *nce, integer 
	*ice, integer *jce, S_fp force)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal a[6000]	/* was [3][2000] */, d__[96000]	/* was [6][
	    2000][8] */, h__;
    static integer j, k, n;
    static doublereal v[6000]	/* was [3][2000] */, x[6000]	/* was [3][
	    2000] */, a0[6000]	/* was [3][2000] */, h2[8];
    static integer j1;
    static doublereal hx2, tmp0, tmp1, tmp2, tol2, vend[6000]	/* was [3][
	    2000] */, xend[6000]	/* was [3][2000] */, vscal[2000], 
	    xscal[2000], errmax;



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



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    --jcen;
    --stat;
    ngf -= 5;
    --rcrit;
    --rphys;
    s -= 4;
    v0 -= 4;
    x0 -= 4;
    --mass;
    --opt;
    --jce;
    --ice;

    /* Function Body */
    tol2 = *tol * *tol;

/* Calculate arrays used to scale the relative error (R^2 for position and */
/* V^2 for velocity). */
    i__1 = *nbod;
    for (k = 2; k <= i__1; ++k) {
	tmp1 = x0[k * 3 + 1] * x0[k * 3 + 1] + x0[k * 3 + 2] * x0[k * 3 + 2] 
		+ x0[k * 3 + 3] * x0[k * 3 + 3];
	tmp2 = v0[k * 3 + 1] * v0[k * 3 + 1] + v0[k * 3 + 2] * v0[k * 3 + 2] 
		+ v0[k * 3 + 3] * v0[k * 3 + 3];
	xscal[k - 1] = 1. / tmp1;
	vscal[k - 1] = 1. / tmp2;
    }

/* Calculate accelerations at the start of the step */
    (*force)(time, &jcen[1], nbod, nbig, &mass[1], &x0[4], &v0[4], &s[4], &
	    rcrit[1], a0, &stat[1], &ngf[5], ngflag, &opt[1], nce, &ice[1], &
	    jce[1]);

L100:

/* For each value of N, do a modified-midpoint integration with 2N substeps */
    for (n = 1; n <= 8; ++n) {
	h__ = *h0 / ((real) n * 2.);
	h2[n - 1] = .25 / (n * n);
	hx2 = h__ * 2.;

	i__1 = *nbod;
	for (k = 2; k <= i__1; ++k) {
	    x[k * 3 - 3] = x0[k * 3 + 1] + h__ * v0[k * 3 + 1];
	    x[k * 3 - 2] = x0[k * 3 + 2] + h__ * v0[k * 3 + 2];
	    x[k * 3 - 1] = x0[k * 3 + 3] + h__ * v0[k * 3 + 3];
	    v[k * 3 - 3] = v0[k * 3 + 1] + h__ * a0[k * 3 - 3];
	    v[k * 3 - 2] = v0[k * 3 + 2] + h__ * a0[k * 3 - 2];
	    v[k * 3 - 1] = v0[k * 3 + 3] + h__ * a0[k * 3 - 1];
	}
	(*force)(time, &jcen[1], nbod, nbig, &mass[1], x, v, &s[4], &rcrit[1],
		 a, &stat[1], &ngf[5], ngflag, &opt[1], nce, &ice[1], &jce[1])
		;
	i__1 = *nbod;
	for (k = 2; k <= i__1; ++k) {
	    xend[k * 3 - 3] = x0[k * 3 + 1] + hx2 * v[k * 3 - 3];
	    xend[k * 3 - 2] = x0[k * 3 + 2] + hx2 * v[k * 3 - 2];
	    xend[k * 3 - 1] = x0[k * 3 + 3] + hx2 * v[k * 3 - 1];
	    vend[k * 3 - 3] = v0[k * 3 + 1] + hx2 * a[k * 3 - 3];
	    vend[k * 3 - 2] = v0[k * 3 + 2] + hx2 * a[k * 3 - 2];
	    vend[k * 3 - 1] = v0[k * 3 + 3] + hx2 * a[k * 3 - 1];
	}

	i__1 = n;
	for (j = 2; j <= i__1; ++j) {
	    (*force)(time, &jcen[1], nbod, nbig, &mass[1], xend, vend, &s[4], 
		    &rcrit[1], a, &stat[1], &ngf[5], ngflag, &opt[1], nce, &
		    ice[1], &jce[1]);
	    i__2 = *nbod;
	    for (k = 2; k <= i__2; ++k) {
		x[k * 3 - 3] += hx2 * vend[k * 3 - 3];
		x[k * 3 - 2] += hx2 * vend[k * 3 - 2];
		x[k * 3 - 1] += hx2 * vend[k * 3 - 1];
		v[k * 3 - 3] += hx2 * a[k * 3 - 3];
		v[k * 3 - 2] += hx2 * a[k * 3 - 2];
		v[k * 3 - 1] += hx2 * a[k * 3 - 1];
	    }
	    (*force)(time, &jcen[1], nbod, nbig, &mass[1], x, v, &s[4], &
		    rcrit[1], a, &stat[1], &ngf[5], ngflag, &opt[1], nce, &
		    ice[1], &jce[1]);
	    i__2 = *nbod;
	    for (k = 2; k <= i__2; ++k) {
		xend[k * 3 - 3] += hx2 * v[k * 3 - 3];
		xend[k * 3 - 2] += hx2 * v[k * 3 - 2];
		xend[k * 3 - 1] += hx2 * v[k * 3 - 1];
		vend[k * 3 - 3] += hx2 * a[k * 3 - 3];
		vend[k * 3 - 2] += hx2 * a[k * 3 - 2];
		vend[k * 3 - 1] += hx2 * a[k * 3 - 1];
	    }
	}

	(*force)(time, &jcen[1], nbod, nbig, &mass[1], xend, vend, &s[4], &
		rcrit[1], a, &stat[1], &ngf[5], ngflag, &opt[1], nce, &ice[1],
		 &jce[1]);

	i__1 = *nbod;
	for (k = 2; k <= i__1; ++k) {
	    d__[(k + n * 2000) * 6 - 12006] = (xend[k * 3 - 3] + x[k * 3 - 3] 
		    + h__ * vend[k * 3 - 3]) * .5;
	    d__[(k + n * 2000) * 6 - 12005] = (xend[k * 3 - 2] + x[k * 3 - 2] 
		    + h__ * vend[k * 3 - 2]) * .5;
	    d__[(k + n * 2000) * 6 - 12004] = (xend[k * 3 - 1] + x[k * 3 - 1] 
		    + h__ * vend[k * 3 - 1]) * .5;
	    d__[(k + n * 2000) * 6 - 12003] = (vend[k * 3 - 3] + v[k * 3 - 3] 
		    + h__ * a[k * 3 - 3]) * .5;
	    d__[(k + n * 2000) * 6 - 12002] = (vend[k * 3 - 2] + v[k * 3 - 2] 
		    + h__ * a[k * 3 - 2]) * .5;
	    d__[(k + n * 2000) * 6 - 12001] = (vend[k * 3 - 1] + v[k * 3 - 1] 
		    + h__ * a[k * 3 - 1]) * .5;
	}

/* Update the D array, used for polynomial extrapolation */
	for (j = n - 1; j >= 1; --j) {
	    j1 = j + 1;
	    tmp0 = 1. / (h2[j - 1] - h2[n - 1]);
	    tmp1 = tmp0 * h2[j1 - 1];
	    tmp2 = tmp0 * h2[n - 1];
	    i__1 = *nbod;
	    for (k = 2; k <= i__1; ++k) {
		d__[(k + j * 2000) * 6 - 12006] = tmp1 * d__[(k + j1 * 2000) *
			 6 - 12006] - tmp2 * d__[(k + j * 2000) * 6 - 12006];
		d__[(k + j * 2000) * 6 - 12005] = tmp1 * d__[(k + j1 * 2000) *
			 6 - 12005] - tmp2 * d__[(k + j * 2000) * 6 - 12005];
		d__[(k + j * 2000) * 6 - 12004] = tmp1 * d__[(k + j1 * 2000) *
			 6 - 12004] - tmp2 * d__[(k + j * 2000) * 6 - 12004];
		d__[(k + j * 2000) * 6 - 12003] = tmp1 * d__[(k + j1 * 2000) *
			 6 - 12003] - tmp2 * d__[(k + j * 2000) * 6 - 12003];
		d__[(k + j * 2000) * 6 - 12002] = tmp1 * d__[(k + j1 * 2000) *
			 6 - 12002] - tmp2 * d__[(k + j * 2000) * 6 - 12002];
		d__[(k + j * 2000) * 6 - 12001] = tmp1 * d__[(k + j1 * 2000) *
			 6 - 12001] - tmp2 * d__[(k + j * 2000) * 6 - 12001];
	    }
	}

/* After several integrations, test the relative error on extrapolated values */
	if (n > 3) {
	    errmax = 0.;

/* Maximum relative position and velocity errors (last D term added) */
	    i__1 = *nbod;
	    for (k = 2; k <= i__1; ++k) {
/* Computing MAX */
		d__1 = d__[(k + 2000) * 6 - 12006] * d__[(k + 2000) * 6 - 
			12006], d__2 = d__[(k + 2000) * 6 - 12005] * d__[(k + 
			2000) * 6 - 12005], d__1 = max(d__1,d__2), d__2 = d__[
			(k + 2000) * 6 - 12004] * d__[(k + 2000) * 6 - 12004];
		tmp1 = max(d__1,d__2);
/* Computing MAX */
		d__1 = d__[(k + 2000) * 6 - 12003] * d__[(k + 2000) * 6 - 
			12003], d__2 = d__[(k + 2000) * 6 - 12002] * d__[(k + 
			2000) * 6 - 12002], d__1 = max(d__1,d__2), d__2 = d__[
			(k + 2000) * 6 - 12001] * d__[(k + 2000) * 6 - 12001];
		tmp2 = max(d__1,d__2);
/* Computing MAX */
		d__1 = errmax, d__2 = tmp1 * xscal[k - 1], d__1 = max(d__1,
			d__2), d__2 = tmp2 * vscal[k - 1];
		errmax = max(d__1,d__2);
	    }

/* If error is smaller than TOL, update position and velocity arrays, and exit */
	    if (errmax <= tol2) {
		i__1 = *nbod;
		for (k = 2; k <= i__1; ++k) {
		    x0[k * 3 + 1] = d__[(k + 2000) * 6 - 12006];
		    x0[k * 3 + 2] = d__[(k + 2000) * 6 - 12005];
		    x0[k * 3 + 3] = d__[(k + 2000) * 6 - 12004];
		    v0[k * 3 + 1] = d__[(k + 2000) * 6 - 12003];
		    v0[k * 3 + 2] = d__[(k + 2000) * 6 - 12002];
		    v0[k * 3 + 3] = d__[(k + 2000) * 6 - 12001];
		}

		i__1 = n;
		for (j = 2; j <= i__1; ++j) {
		    i__2 = *nbod;
		    for (k = 2; k <= i__2; ++k) {
			x0[k * 3 + 1] += d__[(k + j * 2000) * 6 - 12006];
			x0[k * 3 + 2] += d__[(k + j * 2000) * 6 - 12005];
			x0[k * 3 + 3] += d__[(k + j * 2000) * 6 - 12004];
			v0[k * 3 + 1] += d__[(k + j * 2000) * 6 - 12003];
			v0[k * 3 + 2] += d__[(k + j * 2000) * 6 - 12002];
			v0[k * 3 + 3] += d__[(k + j * 2000) * 6 - 12001];
		    }
		}

/* Save the actual stepsize used */
		*hdid = *h0;

/* Recommend a new stepsize for the next call to this subroutine */
		if (n == 8) {
		    *h0 *= .55;
		}
		if (n < 7) {
		    *h0 *= 1.3;
		}
		return 0;
	    }
	}

    }

/* If errors were too large, redo the step with half the previous step size. */
    *h0 *= .5;
    goto L100;

/* ------------------------------------------------------------------------------ */

} /* mdt_bs1__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MDT_BS2.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Integrates NBOD bodies (of which NBIG are Big) for one timestep H0 */
/* using the Bulirsch-Stoer method. The accelerations are calculated using the */
/* subroutine FORCE. The accuracy of the step is approximately determined */
/* by the tolerance parameter TOL. */

/* N.B. This version only works for conservative systems (i.e. force is a */
/* ===  function of position only) !!!! Hence, non-gravitational forces */
/*      and post-Newtonian corrections cannot be used. */

/* N.B. Input/output must be in coordinates with respect to the central body. */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mdt_bs2__(doublereal *time, doublereal *h0, doublereal *
	hdid, doublereal *tol, doublereal *jcen, integer *nbod, integer *nbig,
	 doublereal *mass, doublereal *x0, doublereal *v0, doublereal *s, 
	doublereal *rphys, doublereal *rcrit, doublereal *ngf, integer *stat, 
	integer *dtflag, integer *ngflag, integer *opt, integer *nce, integer 
	*ice, integer *jce, S_fp force)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal a[6000]	/* was [3][2000] */, b[6000]	/* was [3][
	    2000] */, c__[6000]	/* was [3][2000] */, d__[144000]	/* 
	    was [6][2000][12] */, h__;
    static integer j, k, n;
    static doublereal a0[6000]	/* was [3][2000] */, h2[12];
    static integer j1;
    static doublereal hby2, tmp0, tmp1, tmp2, tol2, h2by2, xend[6000]	/* 
	    was [3][2000] */, vscal[2000], xscal[2000], errmax;



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



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    --jcen;
    --stat;
    ngf -= 5;
    --rcrit;
    --rphys;
    s -= 4;
    v0 -= 4;
    x0 -= 4;
    --mass;
    --opt;
    --jce;
    --ice;

    /* Function Body */
    tol2 = *tol * *tol;

/* Calculate arrays used to scale the relative error (R^2 for position and */
/* V^2 for velocity). */
    i__1 = *nbod;
    for (k = 2; k <= i__1; ++k) {
	tmp1 = x0[k * 3 + 1] * x0[k * 3 + 1] + x0[k * 3 + 2] * x0[k * 3 + 2] 
		+ x0[k * 3 + 3] * x0[k * 3 + 3];
	tmp2 = v0[k * 3 + 1] * v0[k * 3 + 1] + v0[k * 3 + 2] * v0[k * 3 + 2] 
		+ v0[k * 3 + 3] * v0[k * 3 + 3];
	xscal[k - 1] = 1. / tmp1;
	vscal[k - 1] = 1. / tmp2;
    }

/* Calculate accelerations at the start of the step */
    (*force)(time, &jcen[1], nbod, nbig, &mass[1], &x0[4], &v0[4], &s[4], &
	    rcrit[1], a0, &stat[1], &ngf[5], ngflag, &opt[1], nce, &ice[1], &
	    jce[1]);

L100:

/* For each value of N, do a modified-midpoint integration with N substeps */
    for (n = 1; n <= 12; ++n) {
	h__ = *h0 / (doublereal) n;
	hby2 = h__ * .5;
	h2[n - 1] = h__ * h__;
	h2by2 = h2[n - 1] * .5;

	i__1 = *nbod;
	for (k = 2; k <= i__1; ++k) {
	    b[k * 3 - 3] = a0[k * 3 - 3] * .5;
	    b[k * 3 - 2] = a0[k * 3 - 2] * .5;
	    b[k * 3 - 1] = a0[k * 3 - 1] * .5;
	    c__[k * 3 - 3] = 0.;
	    c__[k * 3 - 2] = 0.;
	    c__[k * 3 - 1] = 0.;
	    xend[k * 3 - 3] = h2by2 * a0[k * 3 - 3] + h__ * v0[k * 3 + 1] + 
		    x0[k * 3 + 1];
	    xend[k * 3 - 2] = h2by2 * a0[k * 3 - 2] + h__ * v0[k * 3 + 2] + 
		    x0[k * 3 + 2];
	    xend[k * 3 - 1] = h2by2 * a0[k * 3 - 1] + h__ * v0[k * 3 + 3] + 
		    x0[k * 3 + 3];
	}

	i__1 = n;
	for (j = 2; j <= i__1; ++j) {
	    (*force)(time, &jcen[1], nbod, nbig, &mass[1], xend, &v0[4], &s[4]
		    , &rcrit[1], a, &stat[1], &ngf[5], ngflag, &opt[1], nce, &
		    ice[1], &jce[1]);
	    tmp0 = h__ * (doublereal) j;
	    i__2 = *nbod;
	    for (k = 2; k <= i__2; ++k) {
		b[k * 3 - 3] += a[k * 3 - 3];
		b[k * 3 - 2] += a[k * 3 - 2];
		b[k * 3 - 1] += a[k * 3 - 1];
		c__[k * 3 - 3] += b[k * 3 - 3];
		c__[k * 3 - 2] += b[k * 3 - 2];
		c__[k * 3 - 1] += b[k * 3 - 1];
		xend[k * 3 - 3] = h2[n - 1] * c__[k * 3 - 3] + h2by2 * a0[k * 
			3 - 3] + tmp0 * v0[k * 3 + 1] + x0[k * 3 + 1];
		xend[k * 3 - 2] = h2[n - 1] * c__[k * 3 - 2] + h2by2 * a0[k * 
			3 - 2] + tmp0 * v0[k * 3 + 2] + x0[k * 3 + 2];
		xend[k * 3 - 1] = h2[n - 1] * c__[k * 3 - 1] + h2by2 * a0[k * 
			3 - 1] + tmp0 * v0[k * 3 + 3] + x0[k * 3 + 3];
	    }
	}

	(*force)(time, &jcen[1], nbod, nbig, &mass[1], xend, &v0[4], &s[4], &
		rcrit[1], a, &stat[1], &ngf[5], ngflag, &opt[1], nce, &ice[1],
		 &jce[1]);

	i__1 = *nbod;
	for (k = 2; k <= i__1; ++k) {
	    d__[(k + n * 2000) * 6 - 12006] = xend[k * 3 - 3];
	    d__[(k + n * 2000) * 6 - 12005] = xend[k * 3 - 2];
	    d__[(k + n * 2000) * 6 - 12004] = xend[k * 3 - 1];
	    d__[(k + n * 2000) * 6 - 12003] = h__ * b[k * 3 - 3] + hby2 * a[k 
		    * 3 - 3] + v0[k * 3 + 1];
	    d__[(k + n * 2000) * 6 - 12002] = h__ * b[k * 3 - 2] + hby2 * a[k 
		    * 3 - 2] + v0[k * 3 + 2];
	    d__[(k + n * 2000) * 6 - 12001] = h__ * b[k * 3 - 1] + hby2 * a[k 
		    * 3 - 1] + v0[k * 3 + 3];
	}

/* Update the D array, used for polynomial extrapolation */
	for (j = n - 1; j >= 1; --j) {
	    j1 = j + 1;
	    tmp0 = 1. / (h2[j - 1] - h2[n - 1]);
	    tmp1 = tmp0 * h2[j1 - 1];
	    tmp2 = tmp0 * h2[n - 1];
	    i__1 = *nbod;
	    for (k = 2; k <= i__1; ++k) {
		d__[(k + j * 2000) * 6 - 12006] = tmp1 * d__[(k + j1 * 2000) *
			 6 - 12006] - tmp2 * d__[(k + j * 2000) * 6 - 12006];
		d__[(k + j * 2000) * 6 - 12005] = tmp1 * d__[(k + j1 * 2000) *
			 6 - 12005] - tmp2 * d__[(k + j * 2000) * 6 - 12005];
		d__[(k + j * 2000) * 6 - 12004] = tmp1 * d__[(k + j1 * 2000) *
			 6 - 12004] - tmp2 * d__[(k + j * 2000) * 6 - 12004];
		d__[(k + j * 2000) * 6 - 12003] = tmp1 * d__[(k + j1 * 2000) *
			 6 - 12003] - tmp2 * d__[(k + j * 2000) * 6 - 12003];
		d__[(k + j * 2000) * 6 - 12002] = tmp1 * d__[(k + j1 * 2000) *
			 6 - 12002] - tmp2 * d__[(k + j * 2000) * 6 - 12002];
		d__[(k + j * 2000) * 6 - 12001] = tmp1 * d__[(k + j1 * 2000) *
			 6 - 12001] - tmp2 * d__[(k + j * 2000) * 6 - 12001];
	    }
	}

/* After several integrations, test the relative error on extrapolated values */
	if (n > 3) {
	    errmax = 0.;

/* Maximum relative position and velocity errors (last D term added) */
	    i__1 = *nbod;
	    for (k = 2; k <= i__1; ++k) {
/* Computing MAX */
		d__1 = d__[(k + 2000) * 6 - 12006] * d__[(k + 2000) * 6 - 
			12006], d__2 = d__[(k + 2000) * 6 - 12005] * d__[(k + 
			2000) * 6 - 12005], d__1 = max(d__1,d__2), d__2 = d__[
			(k + 2000) * 6 - 12004] * d__[(k + 2000) * 6 - 12004];
		tmp1 = max(d__1,d__2);
/* Computing MAX */
		d__1 = d__[(k + 2000) * 6 - 12003] * d__[(k + 2000) * 6 - 
			12003], d__2 = d__[(k + 2000) * 6 - 12002] * d__[(k + 
			2000) * 6 - 12005], d__1 = max(d__1,d__2), d__2 = d__[
			(k + 2000) * 6 - 12001] * d__[(k + 2000) * 6 - 12001];
		tmp2 = max(d__1,d__2);
/* Computing MAX */
		d__1 = errmax, d__2 = tmp1 * xscal[k - 1], d__1 = max(d__1,
			d__2), d__2 = tmp2 * vscal[k - 1];
		errmax = max(d__1,d__2);
	    }

/* If error is smaller than TOL, update position and velocity arrays and exit */
	    if (errmax <= tol2) {
		i__1 = *nbod;
		for (k = 2; k <= i__1; ++k) {
		    x0[k * 3 + 1] = d__[(k + 2000) * 6 - 12006];
		    x0[k * 3 + 2] = d__[(k + 2000) * 6 - 12005];
		    x0[k * 3 + 3] = d__[(k + 2000) * 6 - 12004];
		    v0[k * 3 + 1] = d__[(k + 2000) * 6 - 12003];
		    v0[k * 3 + 2] = d__[(k + 2000) * 6 - 12002];
		    v0[k * 3 + 3] = d__[(k + 2000) * 6 - 12001];
		}

		i__1 = n;
		for (j = 2; j <= i__1; ++j) {
		    i__2 = *nbod;
		    for (k = 2; k <= i__2; ++k) {
			x0[k * 3 + 1] += d__[(k + j * 2000) * 6 - 12006];
			x0[k * 3 + 2] += d__[(k + j * 2000) * 6 - 12005];
			x0[k * 3 + 3] += d__[(k + j * 2000) * 6 - 12004];
			v0[k * 3 + 1] += d__[(k + j * 2000) * 6 - 12003];
			v0[k * 3 + 2] += d__[(k + j * 2000) * 6 - 12002];
			v0[k * 3 + 3] += d__[(k + j * 2000) * 6 - 12001];
		    }
		}

/* Save the actual stepsize used */
		*hdid = *h0;

/* Recommend a new stepsize for the next call to this subroutine */
		if (n >= 8) {
		    *h0 *= .55;
		}
		if (n < 7) {
		    *h0 *= 1.3;
		}
		return 0;
	    }
	}

    }

/* If errors were too large, redo the step with half the previous step size. */
    *h0 *= .5;
    goto L100;

/* ------------------------------------------------------------------------------ */

} /* mdt_bs2__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MDT_HY.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Integrates NBOD bodies (of which NBIG are Big) for one timestep H */
/* using a second-order hybrid-symplectic integrator algorithm */

/* DTFLAG = 0 implies first ever call to this subroutine, */
/*        = 1 implies first call since number/masses of objects changed. */
/*        = 2 normal call */

/* N.B. Input/output must be in democratic heliocentric coordinates. */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mdt_hy__(doublereal *time, doublereal *tstart, 
	doublereal *h0, doublereal *tol, doublereal *rmax, doublereal *en, 
	doublereal *am, doublereal *jcen, doublereal *rcen, integer *nbod, 
	integer *nbig, doublereal *m, doublereal *x, doublereal *v, 
	doublereal *s, doublereal *rphys, doublereal *rcrit, doublereal *rce, 
	integer *stat, char *id, doublereal *ngf, integer *algor, integer *
	opt, integer *dtflag, integer *ngflag, integer *opflag, integer *
	colflag, integer *nclo, integer *iclo, integer *jclo, doublereal *
	dclo, doublereal *tclo, doublereal *ixvclo, doublereal *jxvclo, char *
	outfile, char *mem, integer *lmem, ftnlen id_len, ftnlen outfile_len, 
	ftnlen mem_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int mfo_hkce__();
    extern /* Subroutine */ int mco_iden__(doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), mdt_hkce__(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, char *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, char *, char *, 
	    integer *, U_fp, ftnlen, ftnlen, ftnlen), mce_snif__(doublereal *,
	     integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *), mfo_user__(doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal a[6000]	/* was [3][2000] */;
    static integer j;
    extern /* Subroutine */ int drift_one__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal v0[6000]	/* was [3][2000] */, x0[6000]	/* was [3][
	    2000] */;
    static integer ce[2000], ice[2000], jce[2000], nce;
    static doublereal hby2, angf[6000]	/* was [3][2000] */, hrec, temp, ausr[
	    6000]	/* was [3][2000] */;
    static integer iflag;
    static doublereal mvsum[3];
    extern /* Subroutine */ int mfo_hy__(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , mfo_ngf__(integer *, doublereal *, doublereal *, doublereal *, 
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

    /* Parameter adjustments */
    --en;
    --am;
    --jcen;
    ngf -= 5;
    id -= 8;
    --stat;
    --rce;
    --rcrit;
    --rphys;
    s -= 4;
    v -= 4;
    x -= 4;
    --m;
    --opt;
    --iclo;
    --jclo;
    --dclo;
    --tclo;
    ixvclo -= 7;
    jxvclo -= 7;
    outfile -= 80;
    mem -= 80;
    --lmem;

    /* Function Body */
    hby2 = *h0 * .5;
    *nclo = 0;
    *colflag = 0;

/* If accelerations from previous call are not valid, calculate them now */
    if (*dtflag != 2) {
	if (*dtflag == 0) {
	    hrec = *h0;
	}
	mfo_hy__(&jcen[1], nbod, nbig, &m[1], &x[4], &rcrit[1], a, &stat[1]);
	*dtflag = 2;
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    angf[j * 3 - 3] = 0.;
	    angf[j * 3 - 2] = 0.;
	    angf[j * 3 - 1] = 0.;
	    ausr[j * 3 - 3] = 0.;
	    ausr[j * 3 - 2] = 0.;
	    ausr[j * 3 - 1] = 0.;
	}
/* If required, apply non-gravitational and user-defined forces */
	if (opt[8] == 1) {
	    mfo_user__(time, &jcen[1], nbod, nbig, &m[1], &x[4], &v[4], ausr);
	}
	if (*ngflag == 1 || *ngflag == 3) {
	    mfo_ngf__(nbod, &x[4], &v[4], angf, &ngf[5]);
	}
    }

/* Advance interaction Hamiltonian for H/2 */
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	v[j * 3 + 1] += hby2 * (angf[j * 3 - 3] + ausr[j * 3 - 3] + a[j * 3 - 
		3]);
	v[j * 3 + 2] += hby2 * (angf[j * 3 - 2] + ausr[j * 3 - 2] + a[j * 3 - 
		2]);
	v[j * 3 + 3] += hby2 * (angf[j * 3 - 1] + ausr[j * 3 - 1] + a[j * 3 - 
		1]);
    }

/* Advance solar Hamiltonian for H/2 */
    mvsum[0] = 0.;
    mvsum[1] = 0.;
    mvsum[2] = 0.;
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	mvsum[0] += m[j] * v[j * 3 + 1];
	mvsum[1] += m[j] * v[j * 3 + 2];
	mvsum[2] += m[j] * v[j * 3 + 3];
    }

    temp = hby2 / m[1];
    mvsum[0] = temp * mvsum[0];
    mvsum[1] = temp * mvsum[1];
    mvsum[2] = temp * mvsum[2];
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	x[j * 3 + 1] += mvsum[0];
	x[j * 3 + 2] += mvsum[1];
	x[j * 3 + 3] += mvsum[2];
    }

/* Save the current coordinates and velocities */
    mco_iden__(time, &jcen[1], nbod, nbig, h0, &m[1], &x[4], &v[4], x0, v0, &
	    ngf[5], ngflag, &opt[1]);

/* Advance H_K for H */
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	iflag = 0;
	drift_one__(&m[1], &x[j * 3 + 1], &x[j * 3 + 2], &x[j * 3 + 3], &v[j *
		 3 + 1], &v[j * 3 + 2], &v[j * 3 + 3], h0, &iflag);
    }

/* Check whether any object separations were < R_CRIT whilst advancing H_K */
    mce_snif__(h0, &c__2, nbod, nbig, x0, v0, &x[4], &v[4], &rcrit[1], ce, &
	    nce, ice, jce);

/* If objects had close encounters, advance H_K using Bulirsch-Stoer instead */
    if (nce > 0) {
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    if (ce[j - 1] != 0) {
		x[j * 3 + 1] = x0[j * 3 - 3];
		x[j * 3 + 2] = x0[j * 3 - 2];
		x[j * 3 + 3] = x0[j * 3 - 1];
		v[j * 3 + 1] = v0[j * 3 - 3];
		v[j * 3 + 2] = v0[j * 3 - 2];
		v[j * 3 + 3] = v0[j * 3 - 1];
	    }
	}
	mdt_hkce__(time, tstart, h0, &hrec, tol, rmax, &en[3], &jcen[1], rcen,
		 nbod, nbig, &m[1], &x[4], &v[4], &s[4], &rphys[1], &rcrit[1],
		 &rce[1], &stat[1], id + 8, &ngf[5], algor, &opt[1], ngflag, 
		colflag, ce, &nce, ice, jce, nclo, &iclo[1], &jclo[1], &dclo[
		1], &tclo[1], &ixvclo[7], &jxvclo[7], outfile + 80, mem + 80, 
		&lmem[1], (U_fp)mfo_hkce__, (ftnlen)8, (ftnlen)80, (ftnlen)80)
		;
    }

/* Advance solar Hamiltonian for H/2 */
    mvsum[0] = 0.;
    mvsum[1] = 0.;
    mvsum[2] = 0.;
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	mvsum[0] += m[j] * v[j * 3 + 1];
	mvsum[1] += m[j] * v[j * 3 + 2];
	mvsum[2] += m[j] * v[j * 3 + 3];
    }

    temp = hby2 / m[1];
    mvsum[0] = temp * mvsum[0];
    mvsum[1] = temp * mvsum[1];
    mvsum[2] = temp * mvsum[2];
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	x[j * 3 + 1] += mvsum[0];
	x[j * 3 + 2] += mvsum[1];
	x[j * 3 + 3] += mvsum[2];
    }

/* Advance interaction Hamiltonian for H/2 */
    mfo_hy__(&jcen[1], nbod, nbig, &m[1], &x[4], &rcrit[1], a, &stat[1]);
    if (opt[8] == 1) {
	mfo_user__(time, &jcen[1], nbod, nbig, &m[1], &x[4], &v[4], ausr);
    }
    if (*ngflag == 1 || *ngflag == 3) {
	mfo_ngf__(nbod, &x[4], &v[4], angf, &ngf[5]);
    }

    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	v[j * 3 + 1] += hby2 * (angf[j * 3 - 3] + ausr[j * 3 - 3] + a[j * 3 - 
		3]);
	v[j * 3 + 2] += hby2 * (angf[j * 3 - 2] + ausr[j * 3 - 2] + a[j * 3 - 
		2]);
	v[j * 3 + 3] += hby2 * (angf[j * 3 - 1] + ausr[j * 3 - 1] + a[j * 3 - 
		1]);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mdt_hy__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MDT_HKCE.FOR    (ErikSoft   1 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Integrates NBOD bodies (of which NBIG are Big) for one timestep H under */
/* the Hamiltonian H_K, including close-encounter terms. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mdt_hkce__(doublereal *time, doublereal *tstart, 
	doublereal *h0, doublereal *hrec, doublereal *tol, doublereal *rmax, 
	doublereal *elost, doublereal *jcen, doublereal *rcen, integer *nbod, 
	integer *nbig, doublereal *m, doublereal *x, doublereal *v, 
	doublereal *s, doublereal *rphy, doublereal *rcrit, doublereal *rce, 
	integer *stat, char *id, doublereal *ngf, integer *algor, integer *
	opt, integer *ngflag, integer *colflag, integer *ce, integer *nce, 
	integer *ice, integer *jce, integer *nclo, integer *iclo, integer *
	jclo, doublereal *dclo, doublereal *tclo, doublereal *ixvclo, 
	doublereal *jxvclo, char *outfile, char *mem, integer *lmem, S_fp 
	force, ftnlen id_len, ftnlen outfile_len, ftnlen mem_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    extern /* Subroutine */ int mco_iden__(doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), mce_coll__(doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, char *, integer *, char *, integer *, 
	    char *, ftnlen, ftnlen, ftnlen);
    static integer nclo_old__;
    extern /* Subroutine */ int mce_stat__(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, char *, char *, integer *, ftnlen, ftnlen);
    static integer i__, j, k;
    static doublereal v0[6000]	/* was [3][2000] */, x0[6000]	/* was [3][
	    2000] */;
    static integer ibs[2000], jbs[2000];
    static doublereal mbs[2000];
    static integer nbs;
    static doublereal sbs[6000]	/* was [3][2000] */, vbs[6000]	/* was [3][
	    2000] */, xbs[6000]	/* was [3][2000] */, tmp0, hdid;
    static char idbs[8*2000];
    static integer ihit[50], jhit[50], chit[50], nhit;
    static doublereal dhit[50], temp, thit[50], thit1;
    static integer iback[2000];
    static doublereal rcebs[2000], ngfbs[8000]	/* was [4][2000] */;
    static integer index[2000], dtflag;
    static doublereal hlocal;
    static integer nbsbig;
    static doublereal tlocal;
    static integer statbs[2000];
    static doublereal rphybs[2000];
    extern /* Subroutine */ int mdt_bs2__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, S_fp);
    static integer nowflag;
    static doublereal rcritbs[2000];



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

/* N.B. Don't set nclo to zero!! */
    /* Parameter adjustments */
    --jcen;
    --ce;
    ngf -= 5;
    id -= 8;
    --stat;
    --rce;
    --rcrit;
    --rphy;
    s -= 4;
    v -= 4;
    x -= 4;
    --m;
    --opt;
    --jce;
    --ice;
    --iclo;
    --jclo;
    --dclo;
    --tclo;
    ixvclo -= 7;
    jxvclo -= 7;
    outfile -= 80;
    mem -= 80;
    --lmem;

    /* Function Body */
    nbs = 1;
    nbsbig = 0;
    mbs[0] = m[1];
    if (*algor == 11) {
	mbs[0] = m[1] + m[2];
    }
    sbs[0] = s[4];
    sbs[1] = s[5];
    sbs[2] = s[6];

/* Put data for close-encounter bodies into local arrays for use with BS routine */
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	if (ce[j] != 0) {
	    ++nbs;
	    if (j <= *nbig) {
		nbsbig = nbs;
	    }
	    mbs[nbs - 1] = m[j];
	    xbs[nbs * 3 - 3] = x[j * 3 + 1];
	    xbs[nbs * 3 - 2] = x[j * 3 + 2];
	    xbs[nbs * 3 - 1] = x[j * 3 + 3];
	    vbs[nbs * 3 - 3] = v[j * 3 + 1];
	    vbs[nbs * 3 - 2] = v[j * 3 + 2];
	    vbs[nbs * 3 - 1] = v[j * 3 + 3];
	    sbs[nbs * 3 - 3] = s[j * 3 + 1];
	    sbs[nbs * 3 - 2] = s[j * 3 + 2];
	    sbs[nbs * 3 - 1] = s[j * 3 + 3];
	    rcebs[nbs - 1] = rce[j];
	    rphybs[nbs - 1] = rphy[j];
	    statbs[nbs - 1] = stat[j];
	    rcritbs[nbs - 1] = rcrit[j];
	    s_copy(idbs + (nbs - 1 << 3), id + (j << 3), (ftnlen)8, (ftnlen)8)
		    ;
	    index[nbs - 1] = j;
	    iback[j - 1] = nbs;
	}
    }

    i__1 = *nce;
    for (k = 1; k <= i__1; ++k) {
	ibs[k - 1] = iback[ice[k] - 1];
	jbs[k - 1] = iback[jce[k] - 1];
    }

    tlocal = 0.;
    hlocal = d_sign(hrec, h0);

/* Begin the Bulirsch-Stoer integration */
L50:
    tmp0 = abs(*h0) - abs(tlocal);
    *hrec = hlocal;
    if (abs(hlocal) > tmp0) {
	hlocal = d_sign(&tmp0, h0);
    }

/* Save old coordinates and integrate */
    mco_iden__(time, &jcen[1], &nbs, &c__0, h0, mbs, xbs, vbs, x0, v0, &ngf[5]
	    , ngflag, &opt[1]);
    mdt_bs2__(time, &hlocal, &hdid, tol, &jcen[1], &nbs, &nbsbig, mbs, xbs, 
	    vbs, sbs, rphybs, rcritbs, ngfbs, statbs, &dtflag, ngflag, &opt[1]
	    , nce, ibs, jbs, (S_fp)force);
    tlocal += hdid;

/* Check for close-encounter minima */
    nclo_old__ = *nclo;
    temp = *time + tlocal;
    mce_stat__(&temp, &hdid, rcen, &nbs, &nbsbig, mbs, x0, v0, xbs, vbs, 
	    rcebs, rphybs, nclo, &iclo[1], &jclo[1], &dclo[1], &tclo[1], &
	    ixvclo[7], &jxvclo[7], &nhit, ihit, jhit, chit, dhit, thit, &
	    thit1, &nowflag, statbs, outfile + 240, mem + 80, &lmem[1], (
	    ftnlen)80, (ftnlen)80);

/* If collisions occurred, resolve the collision and return a flag */
    if (nhit > 0 && opt[2] != 0) {
	i__1 = nhit;
	for (k = 1; k <= i__1; ++k) {
	    if (chit[k - 1] == 1) {
		i__ = ihit[k - 1];
		j = jhit[k - 1];
		mce_coll__(&thit[k - 1], tstart, elost, &jcen[1], &i__, &j, &
			nbs, &nbsbig, mbs, xbs, vbs, sbs, rphybs, statbs, 
			idbs, &opt[1], mem + 80, &lmem[1], outfile + 240, (
			ftnlen)8, (ftnlen)80, (ftnlen)80);
		++(*colflag);
	    }
	}
    }

/* If necessary, continue integrating objects undergoing close encounters */
    if ((tlocal - *h0) * *h0 < 0.) {
	goto L50;
    }

/* Return data for the close-encounter objects to global arrays */
    i__1 = nbs;
    for (k = 2; k <= i__1; ++k) {
	j = index[k - 1];
	m[j] = mbs[k - 1];
	x[j * 3 + 1] = xbs[k * 3 - 3];
	x[j * 3 + 2] = xbs[k * 3 - 2];
	x[j * 3 + 3] = xbs[k * 3 - 1];
	v[j * 3 + 1] = vbs[k * 3 - 3];
	v[j * 3 + 2] = vbs[k * 3 - 2];
	v[j * 3 + 3] = vbs[k * 3 - 1];
	s[j * 3 + 1] = sbs[k * 3 - 3];
	s[j * 3 + 2] = sbs[k * 3 - 2];
	s[j * 3 + 3] = sbs[k * 3 - 1];
	stat[j] = statbs[k - 1];
    }
    i__1 = *nclo;
    for (k = 1; k <= i__1; ++k) {
	iclo[k] = index[iclo[k] - 1];
	jclo[k] = index[jclo[k] - 1];
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mdt_hkce__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MDT_MVS.FOR    (ErikSoft   28 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Integrates NBOD bodies (of which NBIG are Big) for one timestep H */
/* using a second-order mixed-variable symplectic integrator. */

/* DTFLAG = 0 implies first ever call to this subroutine, */
/*        = 1 implies first call since number/masses of objects changed. */
/*        = 2 normal call */

/* N.B. Input/output must be in coordinates with respect to the central body. */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mdt_mvs__(doublereal *time, doublereal *tstart, 
	doublereal *h0, doublereal *tol, doublereal *rmax, doublereal *en, 
	doublereal *am, doublereal *jcen, doublereal *rcen, integer *nbod, 
	integer *nbig, doublereal *m, doublereal *x, doublereal *v, 
	doublereal *s, doublereal *rphys, doublereal *rcrit, doublereal *rce, 
	integer *stat, char *id, doublereal *ngf, integer *algor, integer *
	opt, integer *dtflag, integer *ngflag, integer *opflag, integer *
	colflag, integer *nclo, integer *iclo, integer *jclo, doublereal *
	dclo, doublereal *tclo, doublereal *ixvclo, doublereal *jxvclo, char *
	outfile, char *mem, integer *lmem, ftnlen id_len, ftnlen outfile_len, 
	ftnlen mem_len)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int mco_iden__(doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), mce_stat__(doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    char *, char *, integer *, ftnlen, ftnlen), mfo_user__(doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static doublereal a[6000]	/* was [3][2000] */;
    static integer j;
    extern /* Subroutine */ int drift_one__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal v0[6000]	/* was [3][2000] */, x0[6000]	/* was [3][
	    2000] */, gm[2000], vj[6000]	/* was [3][2000] */, xj[6000]	
	    /* was [3][2000] */, hby2, angf[6000]	/* was [3][2000] */;
    static integer ihit[50], jhit[50], chit[50], nhit;
    static doublereal dhit[50], temp, thit[50], ausr[6000]	/* was [3][
	    2000] */, thit1;
    static integer iflag;
    static doublereal msofar;
    extern /* Subroutine */ int mco_h2j__(doublereal *, doublereal *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), mco_j2h__(doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    mfo_ngf__(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal minside;
    static integer nowflag;
    extern /* Subroutine */ int mfo_mvs__(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;



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
    --en;
    --am;
    --jcen;
    ngf -= 5;
    id -= 8;
    --stat;
    --rce;
    --rcrit;
    --rphys;
    s -= 4;
    v -= 4;
    x -= 4;
    --m;
    --opt;
    --iclo;
    --jclo;
    --dclo;
    --tclo;
    ixvclo -= 7;
    jxvclo -= 7;
    outfile -= 80;
    mem -= 80;
    --lmem;

    /* Function Body */
    hby2 = *h0 * .5;
    *nclo = 0;

/* If accelerations from previous call are not valid, calculate them now, */
/* and also the Jacobi coordinates XJ, and effective central masses GM. */
    if (*dtflag != 2) {
	*dtflag = 2;
	mco_h2j__(time, &jcen[1], nbig, nbig, h0, &m[1], &x[4], &v[4], xj, vj,
		 &ngf[5], ngflag, &opt[1]);
	mfo_mvs__(&jcen[1], nbod, nbig, &m[1], &x[4], xj, a, &stat[1]);

	minside = m[1];
	i__1 = *nbig;
	for (j = 2; j <= i__1; ++j) {
	    msofar = minside + m[j];
	    gm[j - 1] = m[1] * msofar / minside;
	    minside = msofar;
	    angf[j * 3 - 3] = 0.;
	    angf[j * 3 - 2] = 0.;
	    angf[j * 3 - 1] = 0.;
	    ausr[j * 3 - 3] = 0.;
	    ausr[j * 3 - 2] = 0.;
	    ausr[j * 3 - 1] = 0.;
	}
/* If required, apply non-gravitational and user-defined forces */
	if (opt[8] == 1) {
	    mfo_user__(time, &jcen[1], nbod, nbig, &m[1], &x[4], &v[4], ausr);
	}
	if (*ngflag == 1 || *ngflag == 3) {
	    mfo_ngf__(nbod, &x[4], &v[4], angf, &ngf[5]);
	}
    }

/* Advance interaction Hamiltonian for H/2 */
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	v[j * 3 + 1] += hby2 * (angf[j * 3 - 3] + ausr[j * 3 - 3] + a[j * 3 - 
		3]);
	v[j * 3 + 2] += hby2 * (angf[j * 3 - 2] + ausr[j * 3 - 2] + a[j * 3 - 
		2]);
	v[j * 3 + 3] += hby2 * (angf[j * 3 - 1] + ausr[j * 3 - 1] + a[j * 3 - 
		1]);
    }

/* Save current coordinates and velocities */
    mco_iden__(time, &jcen[1], nbod, nbig, h0, &m[1], &x[4], &v[4], x0, v0, &
	    ngf[5], ngflag, &opt[1]);

/* Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies) */
    mco_h2j__(time, &jcen[1], nbig, nbig, h0, &m[1], &x[4], &v[4], xj, vj, &
	    ngf[5], ngflag, &opt[1]);
    i__1 = *nbig;
    for (j = 2; j <= i__1; ++j) {
	iflag = 0;
	drift_one__(&gm[j - 1], &xj[j * 3 - 3], &xj[j * 3 - 2], &xj[j * 3 - 1]
		, &vj[j * 3 - 3], &vj[j * 3 - 2], &vj[j * 3 - 1], h0, &iflag);
    }
    i__1 = *nbod;
    for (j = *nbig + 1; j <= i__1; ++j) {
	iflag = 0;
	drift_one__(&m[1], &x[j * 3 + 1], &x[j * 3 + 2], &x[j * 3 + 3], &v[j *
		 3 + 1], &v[j * 3 + 2], &v[j * 3 + 3], h0, &iflag);
    }
    mco_j2h__(time, &jcen[1], nbig, nbig, h0, &m[1], xj, vj, &x[4], &v[4], &
	    ngf[5], ngflag, &opt[1]);

/* Check for close-encounter minima during drift step */
    temp = *time + *h0;
    mce_stat__(&temp, h0, rcen, nbod, nbig, &m[1], x0, v0, &x[4], &v[4], &rce[
	    1], &rphys[1], nclo, &iclo[1], &jclo[1], &dclo[1], &tclo[1], &
	    ixvclo[7], &jxvclo[7], &nhit, ihit, jhit, chit, dhit, thit, &
	    thit1, &nowflag, &stat[1], outfile + 240, mem + 80, &lmem[1], (
	    ftnlen)80, (ftnlen)80);

/* Advance interaction Hamiltonian for H/2 */
    mfo_mvs__(&jcen[1], nbod, nbig, &m[1], &x[4], xj, a, &stat[1]);
    if (opt[8] == 1) {
	mfo_user__(time, &jcen[1], nbod, nbig, &m[1], &x[4], &v[4], ausr);
    }
    if (*ngflag == 1 || *ngflag == 3) {
	mfo_ngf__(nbod, &x[4], &v[4], angf, &ngf[5]);
    }

    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	v[j * 3 + 1] += hby2 * (angf[j * 3 - 3] + ausr[j * 3 - 3] + a[j * 3 - 
		3]);
	v[j * 3 + 2] += hby2 * (angf[j * 3 - 2] + ausr[j * 3 - 2] + a[j * 3 - 
		2]);
	v[j * 3 + 3] += hby2 * (angf[j * 3 - 1] + ausr[j * 3 - 1] + a[j * 3 - 
		1]);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mdt_mvs__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MDT_RA15.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Integrates NBOD bodies (of which NBIG are Big) for one timestep H0 using */
/* Everhart's RA15 integrator algorithm. The accelerations are calculated */
/* using the subroutine FORCE. The accuracy of the step is approximately */
/* determined by the tolerance parameter TOL. */

/* Based on RADAU by E. Everhart, Physics Department, University of Denver. */
/* Comments giving equation numbers refer to Everhart (1985) ``An Efficient */
/* Integrator that Uses Gauss-Radau Spacings'', in The Dynamics of Comets: */
/* Their Origin and Evolution, p185-202, eds. A. Carusi & G. B. Valsecchi, */
/* pub Reidel. (A listing of the original subroutine is also given in this */
/* paper.) */

/* DTFLAG = 0 implies first ever call to this subroutine, */
/*        = 1 implies first call since number/masses of objects changed. */
/*        = 2 normal call */

/* N.B. Input/output must be in coordinates with respect to the central body. */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mdt_ra15__(doublereal *time, doublereal *t, doublereal *
	tdid, doublereal *tol, doublereal *jcen, integer *nbod, integer *nbig,
	 doublereal *mass, doublereal *x1, doublereal *v1, doublereal *spin, 
	doublereal *rphys, doublereal *rcrit, doublereal *ngf, integer *stat, 
	integer *dtflag, integer *ngflag, integer *opt, integer *nce, integer 
	*ice, integer *jce, S_fp force)
{
    /* Initialized data */

    static doublereal h__[8] = { 0.,.0562625605369221,.1802406917368924,
	    .3526247171131696,.5471536263305554,.7342101772154105,
	    .8853209468390958,.9775206135612875 };
    static doublereal xc[8] = { .5,.1666666666666667,.08333333333333333,.05,
	    .03333333333333333,.02380952380952381,.01785714285714286,
	    .01388888888888889 };
    static doublereal vc[7] = { .5,.3333333333333333,.25,.2,.1666666666666667,
	    .1428571428571429,.125 };

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), d_sign(doublereal *, 
	    doublereal *);

    /* Local variables */
    static doublereal a[6000], b[42000]	/* was [7][6000] */, c__[21], d__[21],
	     e[42000]	/* was [7][6000] */, g[42000]	/* was [7][6000] */;
    static integer j, k, n;
    static doublereal q, r__[28], s[9], v[6000], x[6000], a1[6000], q2, q3, 
	    q4, q5, q6, q7, gk;
    static integer nv;
    static doublereal temp;
    static integer niter;



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


/* Gauss-Radau spacings for substeps within a sequence, for the 15th order */
/* integrator. The sum of the H values should be 3.733333333333333 */

    /* Parameter adjustments */
    --jcen;
    --stat;
    ngf -= 5;
    --rcrit;
    --rphys;
    --spin;
    --v1;
    --x1;
    --mass;
    --opt;
    --jce;
    --ice;

    /* Function Body */

/* Constant coefficients used in series expansions for X and V */
/*  XC: 1/2,  1/6,  1/12, 1/20, 1/30, 1/42, 1/56, 1/72 */
/*  VC: 1/2,  1/3,  1/4,  1/5,  1/6,  1/7,  1/8 */

/* If this is first call to the subroutine, set values of the constant arrays */
/* (R = R21, R31, R32, R41, R42, R43 in Everhart's paper.) */
    if (*dtflag == 0) {
	n = 0;
	for (j = 2; j <= 8; ++j) {
	    i__1 = j - 1;
	    for (k = 1; k <= i__1; ++k) {
		++n;
		r__[n - 1] = 1. / (h__[j - 1] - h__[k - 1]);
	    }
	}

/* Constants to convert between B and G arrays (C = C21, C31, C32, C41, C42...) */
	c__[0] = -h__[1];
	d__[0] = h__[1];
	n = 1;
	for (j = 3; j <= 7; ++j) {
	    ++n;
	    c__[n - 1] = -h__[j - 1] * c__[n - j + 1];
	    d__[n - 1] = h__[1] * d__[n - j + 1];
	    i__1 = j - 1;
	    for (k = 3; k <= i__1; ++k) {
		++n;
		c__[n - 1] = c__[n - j] - h__[j - 1] * c__[n - j + 1];
		d__[n - 1] = d__[n - j] + h__[k - 1] * d__[n - j + 1];
	    }
	    ++n;
	    c__[n - 1] = c__[n - j] - h__[j - 1];
	    d__[n - 1] = d__[n - j] + h__[j - 1];
	}

	*dtflag = 1;
    }

    nv = *nbod * 3;
L100:

/* If this is first call to subroutine since number/masses of objects changed */
/* do 6 iterations and initialize B, E arrays, otherwise do 2 iterations. */
    if (*dtflag == 1) {
	niter = 6;
	i__1 = nv;
	for (j = 4; j <= i__1; ++j) {
	    for (k = 1; k <= 7; ++k) {
		b[k + j * 7 - 8] = 0.;
		e[k + j * 7 - 8] = 0.;
	    }
	}
    } else {
	niter = 2;
    }

/* Calculate forces at the start of the sequence */
    (*force)(time, &jcen[1], nbod, nbig, &mass[1], &x1[1], &v1[1], &spin[1], &
	    rcrit[1], a1, &stat[1], &ngf[5], ngflag, &opt[1], nce, &ice[1], &
	    jce[1]);

/* Find G values from B values predicted at the last call (Eqs. 7 of Everhart) */
    i__1 = nv;
    for (k = 4; k <= i__1; ++k) {
	g[k * 7 - 7] = b[k * 7 - 1] * d__[15] + b[k * 7 - 2] * d__[10] + b[k *
		 7 - 3] * d__[6] + b[k * 7 - 4] * d__[3] + b[k * 7 - 5] * d__[
		1] + b[k * 7 - 6] * d__[0] + b[k * 7 - 7];
	g[k * 7 - 6] = b[k * 7 - 1] * d__[16] + b[k * 7 - 2] * d__[11] + b[k *
		 7 - 3] * d__[7] + b[k * 7 - 4] * d__[4] + b[k * 7 - 5] * d__[
		2] + b[k * 7 - 6];
	g[k * 7 - 5] = b[k * 7 - 1] * d__[17] + b[k * 7 - 2] * d__[12] + b[k *
		 7 - 3] * d__[8] + b[k * 7 - 4] * d__[5] + b[k * 7 - 5];
	g[k * 7 - 4] = b[k * 7 - 1] * d__[18] + b[k * 7 - 2] * d__[13] + b[k *
		 7 - 3] * d__[9] + b[k * 7 - 4];
	g[k * 7 - 3] = b[k * 7 - 1] * d__[19] + b[k * 7 - 2] * d__[14] + b[k *
		 7 - 3];
	g[k * 7 - 2] = b[k * 7 - 1] * d__[20] + b[k * 7 - 2];
	g[k * 7 - 1] = b[k * 7 - 1];
    }

/* ------------------------------------------------------------------------------ */

/*  MAIN  LOOP  STARTS  HERE */

/* For each iteration (six for first call to subroutine, two otherwise)... */
    i__1 = niter;
    for (n = 1; n <= i__1; ++n) {

/* For each substep within a sequence... */
	for (j = 2; j <= 8; ++j) {

/* Calculate position predictors using Eqn. 9 of Everhart */
	    s[0] = *t * h__[j - 1];
	    s[1] = s[0] * s[0] * .5;
	    s[2] = s[1] * h__[j - 1] * .3333333333333333;
	    s[3] = s[2] * h__[j - 1] * .5;
	    s[4] = s[3] * h__[j - 1] * .6;
	    s[5] = s[4] * h__[j - 1] * .6666666666666667;
	    s[6] = s[5] * h__[j - 1] * .7142857142857143;
	    s[7] = s[6] * h__[j - 1] * .75;
	    s[8] = s[7] * h__[j - 1] * .7777777777777778;

	    i__2 = nv;
	    for (k = 4; k <= i__2; ++k) {
		x[k - 1] = s[8] * b[k * 7 - 1] + s[7] * b[k * 7 - 2] + s[6] * 
			b[k * 7 - 3] + s[5] * b[k * 7 - 4] + s[4] * b[k * 7 - 
			5] + s[3] * b[k * 7 - 6] + s[2] * b[k * 7 - 7] + s[1] 
			* a1[k - 1] + s[0] * v1[k] + x1[k];
	    }

/* If necessary, calculate velocity predictors too, from Eqn. 10 of Everhart */
	    if (*ngflag != 0) {
		s[0] = *t * h__[j - 1];
		s[1] = s[0] * h__[j - 1] * .5;
		s[2] = s[1] * h__[j - 1] * .6666666666666667;
		s[3] = s[2] * h__[j - 1] * .75;
		s[4] = s[3] * h__[j - 1] * .8;
		s[5] = s[4] * h__[j - 1] * .8333333333333333;
		s[6] = s[5] * h__[j - 1] * .8571428571428571;
		s[7] = s[6] * h__[j - 1] * .875;

		i__2 = nv;
		for (k = 4; k <= i__2; ++k) {
		    v[k - 1] = s[7] * b[k * 7 - 1] + s[6] * b[k * 7 - 2] + s[
			    5] * b[k * 7 - 3] + s[4] * b[k * 7 - 4] + s[3] * 
			    b[k * 7 - 5] + s[2] * b[k * 7 - 6] + s[1] * b[k * 
			    7 - 7] + s[0] * a1[k - 1] + v1[k];
		}
	    }

/* Calculate forces at the current substep */
	    (*force)(time, &jcen[1], nbod, nbig, &mass[1], x, v, &spin[1], &
		    rcrit[1], a, &stat[1], &ngf[5], ngflag, &opt[1], nce, &
		    ice[1], &jce[1]);

/* Update G values using Eqs. 4 of Everhart, and update B values using Eqs. 5 */
	    if (j == 2) {
		i__2 = nv;
		for (k = 4; k <= i__2; ++k) {
		    temp = g[k * 7 - 7];
		    g[k * 7 - 7] = (a[k - 1] - a1[k - 1]) * r__[0];
		    b[k * 7 - 7] = b[k * 7 - 7] + g[k * 7 - 7] - temp;
		}
		goto L300;
	    }
	    if (j == 3) {
		i__2 = nv;
		for (k = 4; k <= i__2; ++k) {
		    temp = g[k * 7 - 6];
		    gk = a[k - 1] - a1[k - 1];
		    g[k * 7 - 6] = (gk * r__[1] - g[k * 7 - 7]) * r__[2];
		    temp = g[k * 7 - 6] - temp;
		    b[k * 7 - 7] += temp * c__[0];
		    b[k * 7 - 6] += temp;
		}
		goto L300;
	    }
	    if (j == 4) {
		i__2 = nv;
		for (k = 4; k <= i__2; ++k) {
		    temp = g[k * 7 - 5];
		    gk = a[k - 1] - a1[k - 1];
		    g[k * 7 - 5] = ((gk * r__[3] - g[k * 7 - 7]) * r__[4] - g[
			    k * 7 - 6]) * r__[5];
		    temp = g[k * 7 - 5] - temp;
		    b[k * 7 - 7] += temp * c__[1];
		    b[k * 7 - 6] += temp * c__[2];
		    b[k * 7 - 5] += temp;
		}
		goto L300;
	    }
	    if (j == 5) {
		i__2 = nv;
		for (k = 4; k <= i__2; ++k) {
		    temp = g[k * 7 - 4];
		    gk = a[k - 1] - a1[k - 1];
		    g[k * 7 - 4] = (((gk * r__[6] - g[k * 7 - 7]) * r__[7] - 
			    g[k * 7 - 6]) * r__[8] - g[k * 7 - 5]) * r__[9];
		    temp = g[k * 7 - 4] - temp;
		    b[k * 7 - 7] += temp * c__[3];
		    b[k * 7 - 6] += temp * c__[4];
		    b[k * 7 - 5] += temp * c__[5];
		    b[k * 7 - 4] += temp;
		}
		goto L300;
	    }
	    if (j == 6) {
		i__2 = nv;
		for (k = 4; k <= i__2; ++k) {
		    temp = g[k * 7 - 3];
		    gk = a[k - 1] - a1[k - 1];
		    g[k * 7 - 3] = ((((gk * r__[10] - g[k * 7 - 7]) * r__[11] 
			    - g[k * 7 - 6]) * r__[12] - g[k * 7 - 5]) * r__[
			    13] - g[k * 7 - 4]) * r__[14];
		    temp = g[k * 7 - 3] - temp;
		    b[k * 7 - 7] += temp * c__[6];
		    b[k * 7 - 6] += temp * c__[7];
		    b[k * 7 - 5] += temp * c__[8];
		    b[k * 7 - 4] += temp * c__[9];
		    b[k * 7 - 3] += temp;
		}
		goto L300;
	    }
	    if (j == 7) {
		i__2 = nv;
		for (k = 4; k <= i__2; ++k) {
		    temp = g[k * 7 - 2];
		    gk = a[k - 1] - a1[k - 1];
		    g[k * 7 - 2] = (((((gk * r__[15] - g[k * 7 - 7]) * r__[16]
			     - g[k * 7 - 6]) * r__[17] - g[k * 7 - 5]) * r__[
			    18] - g[k * 7 - 4]) * r__[19] - g[k * 7 - 3]) * 
			    r__[20];
		    temp = g[k * 7 - 2] - temp;
		    b[k * 7 - 7] += temp * c__[10];
		    b[k * 7 - 6] += temp * c__[11];
		    b[k * 7 - 5] += temp * c__[12];
		    b[k * 7 - 4] += temp * c__[13];
		    b[k * 7 - 3] += temp * c__[14];
		    b[k * 7 - 2] += temp;
		}
		goto L300;
	    }
	    if (j == 8) {
		i__2 = nv;
		for (k = 4; k <= i__2; ++k) {
		    temp = g[k * 7 - 1];
		    gk = a[k - 1] - a1[k - 1];
		    g[k * 7 - 1] = ((((((gk * r__[21] - g[k * 7 - 7]) * r__[
			    22] - g[k * 7 - 6]) * r__[23] - g[k * 7 - 5]) * 
			    r__[24] - g[k * 7 - 4]) * r__[25] - g[k * 7 - 3]) 
			    * r__[26] - g[k * 7 - 2]) * r__[27];
		    temp = g[k * 7 - 1] - temp;
		    b[k * 7 - 7] += temp * c__[15];
		    b[k * 7 - 6] += temp * c__[16];
		    b[k * 7 - 5] += temp * c__[17];
		    b[k * 7 - 4] += temp * c__[18];
		    b[k * 7 - 3] += temp * c__[19];
		    b[k * 7 - 2] += temp * c__[20];
		    b[k * 7 - 1] += temp;
		}
	    }
L300:
	    ;
	}
    }

/* ------------------------------------------------------------------------------ */

/*  END  OF  MAIN  LOOP */

/* Estimate suitable sequence size for the next call to subroutine (Eqs. 15, 16) */
    temp = 0.;
    i__1 = nv;
    for (k = 4; k <= i__1; ++k) {
/* Computing MAX */
	d__2 = temp, d__3 = (d__1 = b[k * 7 - 1], abs(d__1));
	temp = max(d__2,d__3);
    }
/* Computing 7th power */
    d__1 = abs(*t), d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
    temp /= d__2 * (d__1 * d__1) * 72.;
    *tdid = *t;
    if (temp == 0.) {
	*t = *tdid * 1.4;
    } else {
	d__2 = *tol / temp;
	d__1 = pow_dd(&d__2, &c_b167);
	*t = d_sign(&d__1, tdid);
    }

/* If sequence size for the first subroutine call is too big, go back and redo */
/* the sequence using a smaller size. */
    if (*dtflag == 1 && (d__1 = *t / *tdid, abs(d__1)) < 1.) {
	*t *= .8;
	goto L100;
    }

/* If new sequence size is much bigger than the current one, reduce it */
    if ((d__1 = *t / *tdid, abs(d__1)) > 1.4) {
	*t = *tdid * 1.4;
    }

/* Find new position and velocity values at end of the sequence (Eqs. 11, 12) */
    temp = *tdid * *tdid;
    i__1 = nv;
    for (k = 4; k <= i__1; ++k) {
	x1[k] = (xc[7] * b[k * 7 - 1] + xc[6] * b[k * 7 - 2] + xc[5] * b[k * 
		7 - 3] + xc[4] * b[k * 7 - 4] + xc[3] * b[k * 7 - 5] + xc[2] *
		 b[k * 7 - 6] + xc[1] * b[k * 7 - 7] + xc[0] * a1[k - 1]) * 
		temp + v1[k] * *tdid + x1[k];

	v1[k] = (vc[6] * b[k * 7 - 1] + vc[5] * b[k * 7 - 2] + vc[4] * b[k * 
		7 - 3] + vc[3] * b[k * 7 - 4] + vc[2] * b[k * 7 - 5] + vc[1] *
		 b[k * 7 - 6] + vc[0] * b[k * 7 - 7] + a1[k - 1]) * *tdid + 
		v1[k];
    }

/* Predict new B values to use at the start of the next sequence. The predicted */
/* values from the last call are saved as E. The correction, BD, between the */
/* actual and predicted values of B is applied in advance as a correction. */
    q = *t / *tdid;
    q2 = q * q;
    q3 = q * q2;
    q4 = q2 * q2;
    q5 = q2 * q3;
    q6 = q3 * q3;
    q7 = q3 * q4;

    i__1 = nv;
    for (k = 4; k <= i__1; ++k) {
	s[0] = b[k * 7 - 7] - e[k * 7 - 7];
	s[1] = b[k * 7 - 6] - e[k * 7 - 6];
	s[2] = b[k * 7 - 5] - e[k * 7 - 5];
	s[3] = b[k * 7 - 4] - e[k * 7 - 4];
	s[4] = b[k * 7 - 3] - e[k * 7 - 3];
	s[5] = b[k * 7 - 2] - e[k * 7 - 2];
	s[6] = b[k * 7 - 1] - e[k * 7 - 1];

/* Estimate B values for the next sequence (Eqs. 13 of Everhart). */
	e[k * 7 - 7] = q * (b[k * 7 - 1] * 7. + b[k * 7 - 2] * 6. + b[k * 7 - 
		3] * 5. + b[k * 7 - 4] * 4. + b[k * 7 - 5] * 3. + b[k * 7 - 6]
		 * 2. + b[k * 7 - 7]);
	e[k * 7 - 6] = q2 * (b[k * 7 - 1] * 21. + b[k * 7 - 2] * 15. + b[k * 
		7 - 3] * 10. + b[k * 7 - 4] * 6. + b[k * 7 - 5] * 3. + b[k * 
		7 - 6]);
	e[k * 7 - 5] = q3 * (b[k * 7 - 1] * 35. + b[k * 7 - 2] * 20. + b[k * 
		7 - 3] * 10. + b[k * 7 - 4] * 4. + b[k * 7 - 5]);
	e[k * 7 - 4] = q4 * (b[k * 7 - 1] * 35. + b[k * 7 - 2] * 15. + b[k * 
		7 - 3] * 5. + b[k * 7 - 4]);
	e[k * 7 - 3] = q5 * (b[k * 7 - 1] * 21. + b[k * 7 - 2] * 6. + b[k * 7 
		- 3]);
	e[k * 7 - 2] = q6 * (b[k * 7 - 1] * 7. + b[k * 7 - 2]);
	e[k * 7 - 1] = q7 * b[k * 7 - 1];

	b[k * 7 - 7] = e[k * 7 - 7] + s[0];
	b[k * 7 - 6] = e[k * 7 - 6] + s[1];
	b[k * 7 - 5] = e[k * 7 - 5] + s[2];
	b[k * 7 - 4] = e[k * 7 - 4] + s[3];
	b[k * 7 - 3] = e[k * 7 - 3] + s[4];
	b[k * 7 - 2] = e[k * 7 - 2] + s[5];
	b[k * 7 - 1] = e[k * 7 - 1] + s[6];
    }
    *dtflag = 2;

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mdt_ra15__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MFO_ALL.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates accelerations on a set of NBOD bodies (of which NBIG are Big) */
/* due to Newtonian gravitational perturbations, post-Newtonian */
/* corrections (if required), cometary non-gravitational forces (if required) */
/* and user-defined forces (if required). */

/* N.B. Input/output must be in coordinates with respect to the central body. */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mfo_all__(doublereal *time, doublereal *jcen, integer *
	nbod, integer *nbig, doublereal *m, doublereal *x, doublereal *v, 
	doublereal *s, doublereal *rcrit, doublereal *a, integer *stat, 
	doublereal *ngf, integer *ngflag, integer *opt, integer *nce, integer 
	*ice, integer *jce)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int mfo_grav__(integer *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *), mfo_user__(
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer j;
    static doublereal acen[3], acor[6000]	/* was [3][2000] */;
    extern /* Subroutine */ int mfo_pn__(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), mfo_pr__(integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *), mfo_ngf__(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), mfo_obl__(doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);



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

/* Newtonian gravitational forces */
    /* Parameter adjustments */
    --jcen;
    ngf -= 5;
    --stat;
    a -= 4;
    --rcrit;
    s -= 4;
    v -= 4;
    x -= 4;
    --m;
    --opt;
    --jce;
    --ice;

    /* Function Body */
    mfo_grav__(nbod, nbig, &m[1], &x[4], &v[4], &a[4], &stat[1]);

/* Correct for oblateness of the central body */
    if (jcen[1] != 0. || jcen[2] != 0. || jcen[3] != 0.) {
	mfo_obl__(&jcen[1], nbod, &m[1], &x[4], acor, acen);
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    a[j * 3 + 1] += acor[j * 3 - 3] - acen[0];
	    a[j * 3 + 2] += acor[j * 3 - 2] - acen[1];
	    a[j * 3 + 3] += acor[j * 3 - 1] - acen[2];
	}
    }

/* Include non-gravitational (cometary jet) accelerations if necessary */
    if (*ngflag == 1 || *ngflag == 3) {
	mfo_ngf__(nbod, &x[4], &v[4], acor, &ngf[5]);
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    a[j * 3 + 1] += acor[j * 3 - 3];
	    a[j * 3 + 2] += acor[j * 3 - 2];
	    a[j * 3 + 3] += acor[j * 3 - 1];
	}
    }

/* Include radiation pressure/Poynting-Robertson drag if necessary */
    if (*ngflag == 2 || *ngflag == 3) {
	mfo_pr__(nbod, nbig, &m[1], &x[4], &v[4], acor, &ngf[5]);
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    a[j * 3 + 1] += acor[j * 3 - 3];
	    a[j * 3 + 2] += acor[j * 3 - 2];
	    a[j * 3 + 3] += acor[j * 3 - 1];
	}
    }

/* Include post-Newtonian corrections if required */
    if (opt[7] == 1) {
	mfo_pn__(nbod, nbig, &m[1], &x[4], &v[4], acor);
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    a[j * 3 + 1] += acor[j * 3 - 3];
	    a[j * 3 + 2] += acor[j * 3 - 2];
	    a[j * 3 + 3] += acor[j * 3 - 1];
	}
    }

/* Include user-defined accelerations if required */
    if (opt[8] == 1) {
	mfo_user__(time, &jcen[1], nbod, nbig, &m[1], &x[4], &v[4], acor);
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    a[j * 3 + 1] += acor[j * 3 - 3];
	    a[j * 3 + 2] += acor[j * 3 - 2];
	    a[j * 3 + 3] += acor[j * 3 - 1];
	}
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mfo_all__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MFO_GRAV.FOR    (ErikSoft   3 October 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates accelerations on a set of NBOD bodies (NBIG of which are Big) */
/* due to gravitational perturbations by all the other bodies, except that */
/* Small bodies do not interact with one another. */

/* The positions and velocities are stored in arrays X, V with the format */
/* (x,y,z) and (vx,vy,vz) for each object in succession. The accelerations */
/* are stored in the array A (ax,ay,az). */

/* N.B. All coordinates and velocities must be with respect to central body!!!! */
/* === */
/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mfo_grav__(integer *nbod, integer *nbig, doublereal *m, 
	doublereal *x, doublereal *v, doublereal *a, integer *stat)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal s2, r3[2000], dx, dy, dz, sx, sy, sz, s_1__, s_3__, 
	    tmp1, tmp2;



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
    --stat;
    a -= 4;
    v -= 4;
    x -= 4;
    --m;

    /* Function Body */
    sx = 0.;
    sy = 0.;
    sz = 0.;
    i__1 = *nbod;
    for (i__ = 2; i__ <= i__1; ++i__) {
	a[i__ * 3 + 1] = 0.;
	a[i__ * 3 + 2] = 0.;
	a[i__ * 3 + 3] = 0.;
	s2 = x[i__ * 3 + 1] * x[i__ * 3 + 1] + x[i__ * 3 + 2] * x[i__ * 3 + 2]
		 + x[i__ * 3 + 3] * x[i__ * 3 + 3];
	s_1__ = 1. / sqrt(s2);
	r3[i__ - 1] = s_1__ * s_1__ * s_1__;
    }

    i__1 = *nbod;
    for (i__ = 2; i__ <= i__1; ++i__) {
	tmp1 = m[i__] * r3[i__ - 1];
	sx -= tmp1 * x[i__ * 3 + 1];
	sy -= tmp1 * x[i__ * 3 + 2];
	sz -= tmp1 * x[i__ * 3 + 3];
    }

/* Direct terms */
    i__1 = *nbig;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = *nbod;
	for (j = i__ + 1; j <= i__2; ++j) {
	    dx = x[j * 3 + 1] - x[i__ * 3 + 1];
	    dy = x[j * 3 + 2] - x[i__ * 3 + 2];
	    dz = x[j * 3 + 3] - x[i__ * 3 + 3];
	    s2 = dx * dx + dy * dy + dz * dz;
	    s_1__ = 1. / sqrt(s2);
	    s_3__ = s_1__ * s_1__ * s_1__;
	    tmp1 = s_3__ * m[i__];
	    tmp2 = s_3__ * m[j];
	    a[j * 3 + 1] -= tmp1 * dx;
	    a[j * 3 + 2] -= tmp1 * dy;
	    a[j * 3 + 3] -= tmp1 * dz;
	    a[i__ * 3 + 1] += tmp2 * dx;
	    a[i__ * 3 + 2] += tmp2 * dy;
	    a[i__ * 3 + 3] += tmp2 * dz;
	}
    }

/* Indirect terms (add these on last to reduce roundoff error) */
    i__1 = *nbod;
    for (i__ = 2; i__ <= i__1; ++i__) {
	tmp1 = m[1] * r3[i__ - 1];
	a[i__ * 3 + 1] = a[i__ * 3 + 1] + sx - tmp1 * x[i__ * 3 + 1];
	a[i__ * 3 + 2] = a[i__ * 3 + 2] + sy - tmp1 * x[i__ * 3 + 2];
	a[i__ * 3 + 3] = a[i__ * 3 + 3] + sz - tmp1 * x[i__ * 3 + 3];
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mfo_grav__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*     MFO_DRCT.FOR    (ErikSoft   27 February 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates direct accelerations between bodies in the interaction part */
/* of the Hamiltonian of a symplectic integrator that partitions close */
/* encounter terms (e.g. hybrid symplectic algorithms or SyMBA). */
/* The routine calculates accelerations between all pairs of bodies with */
/* indices I >= I0. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mfo_drct__(integer *i0, integer *nbod, integer *nbig, 
	doublereal *m, doublereal *x, doublereal *rcrit, doublereal *a, 
	integer *stat)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal q, s, q2, q3, s2, q4, q5, rc, dx, dy, dz, s_1__, s_3__, 
	    rc2, tmp2, faci, facj;



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
    --stat;
    a -= 4;
    --rcrit;
    x -= 4;
    --m;

    /* Function Body */
    if (*i0 <= 0) {
	*i0 = 2;
    }

    i__1 = *nbig;
    for (i__ = *i0; i__ <= i__1; ++i__) {
	i__2 = *nbod;
	for (j = i__ + 1; j <= i__2; ++j) {
	    dx = x[j * 3 + 1] - x[i__ * 3 + 1];
	    dy = x[j * 3 + 2] - x[i__ * 3 + 2];
	    dz = x[j * 3 + 3] - x[i__ * 3 + 3];
	    s2 = dx * dx + dy * dy + dz * dz;
/* Computing MAX */
	    d__1 = rcrit[i__], d__2 = rcrit[j];
	    rc = max(d__1,d__2);
	    rc2 = rc * rc;

	    if (s2 >= rc2) {
		s_1__ = 1. / sqrt(s2);
		tmp2 = s_1__ * s_1__ * s_1__;
	    } else if (s2 <= rc2 * .01f) {
		tmp2 = 0.;
	    } else {
		s_1__ = 1. / sqrt(s2);
		s = 1. / s_1__;
		s_3__ = s_1__ * s_1__ * s_1__;
		q = (s - rc * .1) / (rc * .9);
		q2 = q * q;
		q3 = q * q2;
		q4 = q2 * q2;
		q5 = q2 * q3;
		tmp2 = (q3 * 10. - q4 * 15. + q5 * 6.) * s_3__;
	    }

	    faci = tmp2 * m[i__];
	    facj = tmp2 * m[j];
	    a[j * 3 + 1] -= faci * dx;
	    a[j * 3 + 2] -= faci * dy;
	    a[j * 3 + 3] -= faci * dz;
	    a[i__ * 3 + 1] += facj * dx;
	    a[i__ * 3 + 2] += facj * dy;
	    a[i__ * 3 + 3] += facj * dz;
	}
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mfo_drct__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*     MFO_HY.FOR    (ErikSoft   2 October 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates accelerations due to the Interaction part of the Hamiltonian */
/* of a hybrid symplectic integrator for a set of NBOD bodies (NBIG of which */
/* are Big), where Small bodies do not interact with one another. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mfo_hy__(doublereal *jcen, integer *nbod, integer *nbig, 
	doublereal *m, doublereal *x, doublereal *rcrit, doublereal *a, 
	integer *stat)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int mfo_drct__(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;
    static integer k;
    static doublereal acen[3], aobl[6000]	/* was [3][2000] */;
    extern /* Subroutine */ int mfo_obl__(doublereal *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *);



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

/* Initialize accelerations to zero */
    /* Parameter adjustments */
    --jcen;
    --stat;
    a -= 4;
    --rcrit;
    x -= 4;
    --m;

    /* Function Body */
    i__1 = *nbod;
    for (k = 1; k <= i__1; ++k) {
	a[k * 3 + 1] = 0.;
	a[k * 3 + 2] = 0.;
	a[k * 3 + 3] = 0.;
    }

/* Calculate direct terms */
    mfo_drct__(&c__2, nbod, nbig, &m[1], &x[4], &rcrit[1], &a[4], &stat[1]);

/* Add accelerations due to oblateness of the central body */
    if (jcen[1] != 0. || jcen[2] != 0. || jcen[3] != 0.) {
	mfo_obl__(&jcen[1], nbod, &m[1], &x[4], aobl, acen);
	i__1 = *nbod;
	for (k = 2; k <= i__1; ++k) {
	    a[k * 3 + 1] = a[k * 3 + 1] + aobl[k * 3 - 3] - acen[0];
	    a[k * 3 + 2] = a[k * 3 + 2] + aobl[k * 3 - 2] - acen[1];
	    a[k * 3 + 3] = a[k * 3 + 3] + aobl[k * 3 - 1] - acen[2];
	}
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mfo_hy__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*     MFO_HKCE.FOR    (ErikSoft   27 February 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates accelerations due to the Keplerian part of the Hamiltonian */
/* of a hybrid symplectic integrator, when close encounters are taking place, */
/* for a set of NBOD bodies (NBIG of which are Big). Note that Small bodies */
/* do not interact with one another. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mfo_hkce__(doublereal *time, doublereal *jcen, integer *
	nbod, integer *nbig, doublereal *m, doublereal *x, doublereal *v, 
	doublereal *spin, doublereal *rcrit, doublereal *a, integer *stat, 
	doublereal *ngf, integer *ngflag, integer *opt, integer *nce, integer 
	*ice, integer *jce)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal q, s, q2, q3, s2, q4, q5, rc, dx, dy, dz, s_1__, s_3__, 
	    rc2, tmp2, faci, facj;



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

/* Initialize accelerations */
    /* Parameter adjustments */
    --jcen;
    ngf -= 5;
    --stat;
    a -= 4;
    --rcrit;
    spin -= 4;
    v -= 4;
    x -= 4;
    --m;
    --opt;
    --jce;
    --ice;

    /* Function Body */
    i__1 = *nbod;
    for (j = 1; j <= i__1; ++j) {
	a[j * 3 + 1] = 0.;
	a[j * 3 + 2] = 0.;
	a[j * 3 + 3] = 0.;
    }

/* Direct terms */
    i__1 = *nce;
    for (k = 1; k <= i__1; ++k) {
	i__ = ice[k];
	j = jce[k];
	dx = x[j * 3 + 1] - x[i__ * 3 + 1];
	dy = x[j * 3 + 2] - x[i__ * 3 + 2];
	dz = x[j * 3 + 3] - x[i__ * 3 + 3];
	s2 = dx * dx + dy * dy + dz * dz;
/* Computing MAX */
	d__1 = rcrit[i__], d__2 = rcrit[j];
	rc = max(d__1,d__2);
	rc2 = rc * rc;

	if (s2 < rc2) {
	    s_1__ = 1. / sqrt(s2);
	    s_3__ = s_1__ * s_1__ * s_1__;
	    if (s2 <= rc2 * .01f) {
		tmp2 = s_3__;
	    } else {
		s = 1. / s_1__;
		q = (s - rc * .1) / (rc * .9);
		q2 = q * q;
		q3 = q * q2;
		q4 = q2 * q2;
		q5 = q2 * q3;
		tmp2 = (1. - q3 * 10. + q4 * 15. - q5 * 6.) * s_3__;
	    }

	    faci = tmp2 * m[i__];
	    facj = tmp2 * m[j];
	    a[j * 3 + 1] -= faci * dx;
	    a[j * 3 + 2] -= faci * dy;
	    a[j * 3 + 3] -= faci * dz;
	    a[i__ * 3 + 1] += facj * dx;
	    a[i__ * 3 + 2] += facj * dy;
	    a[i__ * 3 + 3] += facj * dz;
	}
    }

/* Solar terms */
    i__1 = *nbod;
    for (i__ = 2; i__ <= i__1; ++i__) {
	s2 = x[i__ * 3 + 1] * x[i__ * 3 + 1] + x[i__ * 3 + 2] * x[i__ * 3 + 2]
		 + x[i__ * 3 + 3] * x[i__ * 3 + 3];
	s_1__ = 1. / sqrt(s2);
	tmp2 = m[1] * s_1__ * s_1__ * s_1__;
	a[i__ * 3 + 1] -= tmp2 * x[i__ * 3 + 1];
	a[i__ * 3 + 2] -= tmp2 * x[i__ * 3 + 2];
	a[i__ * 3 + 3] -= tmp2 * x[i__ * 3 + 3];
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mfo_hkce__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*     MFO_MVS.FOR    (ErikSoft   2 October 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates accelerations on a set of NBOD bodies (of which NBIG are Big) */
/* due to gravitational perturbations by all the other bodies. */
/* This routine is designed for use with a mixed-variable symplectic */
/* integrator using Jacobi coordinates. */

/* Based upon routines from Levison and Duncan's SWIFT integrator. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mfo_mvs__(doublereal *jcen, integer *nbod, integer *nbig,
	 doublereal *m, doublereal *x, doublereal *xj, doublereal *a, integer 
	*stat)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal q, r__, a0[3], a1[6000]	/* was [3][2000] */, a2[6000]	
	    /* was [3][2000] */, a3[6000]	/* was [3][2000] */;
    static integer k1;
    static doublereal q2, r2, s2, r3, q3, q4, q5, q6, q7, dx, dy, dz, rj, 
	    s_1__, s_3__, rj2, rj3, fac0, fac1, fac2, a0tp[3], fac12, faci, 
	    facj, acen[3], aobl[6000]	/* was [3][2000] */;
    extern /* Subroutine */ int mfo_obl__(doublereal *, integer *, doublereal 
	    *, doublereal *, doublereal *, doublereal *);
    static doublereal minside;



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

/* Initialize variables */
    /* Parameter adjustments */
    --jcen;
    --stat;
    a -= 4;
    xj -= 4;
    x -= 4;
    --m;

    /* Function Body */
    a0[0] = 0.;
    a0[1] = 0.;
    a0[2] = 0.;
    a1[3] = 0.;
    a1[4] = 0.;
    a1[5] = 0.;
    a2[3] = 0.;
    a2[4] = 0.;
    a2[5] = 0.;
    minside = 0.;

/* Calculate acceleration terms */
    i__1 = *nbig;
    for (k = 3; k <= i__1; ++k) {
	k1 = k - 1;
	minside += m[k1];
	r2 = x[k * 3 + 1] * x[k * 3 + 1] + x[k * 3 + 2] * x[k * 3 + 2] + x[k *
		 3 + 3] * x[k * 3 + 3];
	rj2 = xj[k * 3 + 1] * xj[k * 3 + 1] + xj[k * 3 + 2] * xj[k * 3 + 2] + 
		xj[k * 3 + 3] * xj[k * 3 + 3];
	r__ = 1. / sqrt(r2);
	rj = 1. / sqrt(rj2);
	r3 = r__ * r__ * r__;
	rj3 = rj * rj * rj;

	fac0 = m[k] * r3;
	fac12 = m[1] * rj3;
	fac2 = m[k] * fac12 / (minside + m[1]);
	q = (r2 - rj2) * .5 / rj2;
	q2 = q * q;
	q3 = q * q2;
	q4 = q2 * q2;
	q5 = q2 * q3;
	q6 = q3 * q3;
	q7 = q3 * q4;
	fac1 = q7 * 402.1875 - q6 * 187.6875 + q5 * 86.625 - q4 * 39.375 + q3 
		* 17.5 - q2 * 7.5 + q * 3. - 1.;

/* Add to A0 term */
	a0[0] -= fac0 * x[k * 3 + 1];
	a0[1] -= fac0 * x[k * 3 + 2];
	a0[2] -= fac0 * x[k * 3 + 3];

/* Calculate A1 for this body */
	a1[k * 3 - 3] = fac12 * (xj[k * 3 + 1] + fac1 * x[k * 3 + 1]);
	a1[k * 3 - 2] = fac12 * (xj[k * 3 + 2] + fac1 * x[k * 3 + 2]);
	a1[k * 3 - 1] = fac12 * (xj[k * 3 + 3] + fac1 * x[k * 3 + 3]);

/* Calculate A2 for this body */
	a2[k * 3 - 3] = a2[k1 * 3 - 3] + fac2 * xj[k * 3 + 1];
	a2[k * 3 - 2] = a2[k1 * 3 - 2] + fac2 * xj[k * 3 + 2];
	a2[k * 3 - 1] = a2[k1 * 3 - 1] + fac2 * xj[k * 3 + 3];
    }

    r2 = x[7] * x[7] + x[8] * x[8] + x[9] * x[9];
    r__ = 1. / sqrt(r2);
    r3 = r__ * r__ * r__;
    fac0 = m[2] * r3;
    a0tp[0] = a0[0] - fac0 * x[7];
    a0tp[1] = a0[1] - fac0 * x[8];
    a0tp[2] = a0[2] - fac0 * x[9];

/* Calculate A3 (direct terms) */
    i__1 = *nbod;
    for (k = 2; k <= i__1; ++k) {
	a3[k * 3 - 3] = 0.;
	a3[k * 3 - 2] = 0.;
	a3[k * 3 - 1] = 0.;
    }
    i__1 = *nbig;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = *nbig;
	for (j = i__ + 1; j <= i__2; ++j) {
	    dx = x[j * 3 + 1] - x[i__ * 3 + 1];
	    dy = x[j * 3 + 2] - x[i__ * 3 + 2];
	    dz = x[j * 3 + 3] - x[i__ * 3 + 3];
	    s2 = dx * dx + dy * dy + dz * dz;
	    s_1__ = 1. / sqrt(s2);
	    s_3__ = s_1__ * s_1__ * s_1__;
	    faci = m[i__] * s_3__;
	    facj = m[j] * s_3__;
	    a3[j * 3 - 3] -= faci * dx;
	    a3[j * 3 - 2] -= faci * dy;
	    a3[j * 3 - 1] -= faci * dz;
	    a3[i__ * 3 - 3] += facj * dx;
	    a3[i__ * 3 - 2] += facj * dy;
	    a3[i__ * 3 - 1] += facj * dz;
	}

	i__2 = *nbod;
	for (j = *nbig + 1; j <= i__2; ++j) {
	    dx = x[j * 3 + 1] - x[i__ * 3 + 1];
	    dy = x[j * 3 + 2] - x[i__ * 3 + 2];
	    dz = x[j * 3 + 3] - x[i__ * 3 + 3];
	    s2 = dx * dx + dy * dy + dz * dz;
	    s_1__ = 1. / sqrt(s2);
	    s_3__ = s_1__ * s_1__ * s_1__;
	    faci = m[i__] * s_3__;
	    a3[j * 3 - 3] -= faci * dx;
	    a3[j * 3 - 2] -= faci * dy;
	    a3[j * 3 - 1] -= faci * dz;
	}
    }

/* Big-body accelerations */
    i__1 = *nbig;
    for (k = 2; k <= i__1; ++k) {
	a[k * 3 + 1] = a0[0] + a1[k * 3 - 3] + a2[k * 3 - 3] + a3[k * 3 - 3];
	a[k * 3 + 2] = a0[1] + a1[k * 3 - 2] + a2[k * 3 - 2] + a3[k * 3 - 2];
	a[k * 3 + 3] = a0[2] + a1[k * 3 - 1] + a2[k * 3 - 1] + a3[k * 3 - 1];
    }

/* Small-body accelerations */
    i__1 = *nbod;
    for (k = *nbig + 1; k <= i__1; ++k) {
	a[k * 3 + 1] = a0tp[0] + a3[k * 3 - 3];
	a[k * 3 + 2] = a0tp[1] + a3[k * 3 - 2];
	a[k * 3 + 3] = a0tp[2] + a3[k * 3 - 1];
    }

/* Correct for oblateness of the central body */
    if (jcen[1] != 0. || jcen[2] != 0. || jcen[3] != 0.) {
	mfo_obl__(&jcen[1], nbod, &m[1], &x[4], aobl, acen);
	i__1 = *nbod;
	for (k = 2; k <= i__1; ++k) {
	    a[k * 3 + 1] += aobl[k * 3 - 3] - acen[0];
	    a[k * 3 + 2] += aobl[k * 3 - 2] - acen[1];
	    a[k * 3 + 3] += aobl[k * 3 - 1] - acen[2];
	}
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mfo_mvs__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MFO_NGF.FOR    (ErikSoft  29 November 1999) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates accelerations on a set of NBOD bodies due to cometary */
/* non-gravitational jet forces. The positions and velocities are stored in */
/* arrays X, V with the format (x,y,z) and (vx,vy,vz) for each object in */
/* succession. The accelerations are stored in the array A (ax,ay,az). The */
/* non-gravitational accelerations follow a force law described by Marsden */
/* et al. (1973) Astron. J. 211-225, with magnitude determined by the */
/* parameters NGF(1,2,3) for each object. */

/* N.B. All coordinates and velocities must be with respect to central body!!!! */
/* === */
/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mfo_ngf__(integer *nbod, doublereal *x, doublereal *v, 
	doublereal *a, doublereal *ngf)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal g;
    static integer j;
    static doublereal q, r__, a1, a2, a3, r2, nx, ny, rv, nz, tx, ty, tz;



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
    ngf -= 5;
    a -= 4;
    v -= 4;
    x -= 4;

    /* Function Body */
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	r2 = x[j * 3 + 1] * x[j * 3 + 1] + x[j * 3 + 2] * x[j * 3 + 2] + x[j *
		 3 + 3] * x[j * 3 + 3];

/* Only calculate accelerations if body is close to the Sun (R < 9.36 AU), */
/* or if the non-gravitational force parameters are exceptionally large. */
	if (r2 < 88. || (d__1 = ngf[(j << 2) + 1], abs(d__1)) > 1e-7 || (d__2 
		= ngf[(j << 2) + 2], abs(d__2)) > 1e-7 || (d__3 = ngf[(j << 2)
		 + 3], abs(d__3)) > 1e-7) {
	    r__ = sqrt(r2);
	    rv = x[j * 3 + 1] * v[j * 3 + 1] + x[j * 3 + 2] * v[j * 3 + 2] + 
		    x[j * 3 + 3] * v[j * 3 + 3];

/* Calculate Q = R / R0, where R0 = 2.808 AU */
	    q = r__ * .3561253561253561;
	    d__1 = pow_dd(&q, &c_b178) + 1.;
	    g = pow_dd(&q, &c_b176) * .111262 * pow_dd(&d__1, &c_b177);

/* Within-orbital-plane transverse vector components */
	    tx = r2 * v[j * 3 + 1] - rv * x[j * 3 + 1];
	    ty = r2 * v[j * 3 + 2] - rv * x[j * 3 + 2];
	    tz = r2 * v[j * 3 + 3] - rv * x[j * 3 + 3];

/* Orbit-normal vector components */
	    nx = x[j * 3 + 2] * v[j * 3 + 3] - x[j * 3 + 3] * v[j * 3 + 2];
	    ny = x[j * 3 + 3] * v[j * 3 + 1] - x[j * 3 + 1] * v[j * 3 + 3];
	    nz = x[j * 3 + 1] * v[j * 3 + 2] - x[j * 3 + 2] * v[j * 3 + 1];

/* Multiplication factors */
	    a1 = ngf[(j << 2) + 1] * g / r__;
	    a2 = ngf[(j << 2) + 2] * g / sqrt(tx * tx + ty * ty + tz * tz);
	    a3 = ngf[(j << 2) + 3] * g / sqrt(nx * nx + ny * ny + nz * nz);

/* X,Y and Z components of non-gravitational acceleration */
	    a[j * 3 + 1] = a1 * x[j * 3 + 1] + a2 * tx + a3 * nx;
	    a[j * 3 + 2] = a1 * x[j * 3 + 2] + a2 * ty + a3 * ny;
	    a[j * 3 + 3] = a1 * x[j * 3 + 3] + a2 * tz + a3 * nz;
	} else {
	    a[j * 3 + 1] = 0.;
	    a[j * 3 + 2] = 0.;
	    a[j * 3 + 3] = 0.;
	}
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mfo_ngf__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*     MFO_OBL.FOR    (ErikSoft   2 October 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates barycentric accelerations of NBOD bodies due to oblateness of */
/* the central body. Also returns the corresponding barycentric acceleration */
/* of the central body. */

/* N.B. All coordinates must be with respect to the central body!!!! */
/* === */
/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mfo_obl__(doublereal *jcen, integer *nbod, doublereal *m,
	 doublereal *x, doublereal *a, doublereal *acen)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal r2, u2, u4, u6, r_1__, r_2__, r_3__, jr2, jr4, jr6, 
	    tmp1, tmp2, tmp3, tmp4;



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
    --jcen;
    a -= 4;
    x -= 4;
    --m;
    --acen;

    /* Function Body */
    acen[1] = 0.;
    acen[2] = 0.;
    acen[3] = 0.;

    i__1 = *nbod;
    for (i__ = 2; i__ <= i__1; ++i__) {

/* Calculate barycentric accelerations on the objects */
	r2 = x[i__ * 3 + 1] * x[i__ * 3 + 1] + x[i__ * 3 + 2] * x[i__ * 3 + 2]
		 + x[i__ * 3 + 3] * x[i__ * 3 + 3];
	r_1__ = 1. / sqrt(r2);
	r_2__ = r_1__ * r_1__;
	r_3__ = r_2__ * r_1__;
	jr2 = jcen[1] * r_2__;
	jr4 = jcen[2] * r_2__ * r_2__;
	jr6 = jcen[3] * r_2__ * r_2__ * r_2__;
	u2 = x[i__ * 3 + 3] * x[i__ * 3 + 3] * r_2__;
	u4 = u2 * u2;
	u6 = u4 * u2;

	tmp1 = m[1] * r_3__;
	tmp2 = jr2 * (u2 * 7.5 - 1.5) + jr4 * (u4 * 39.375 - u2 * 26.25 + 
		1.875) + jr6 * (u6 * 187.6875 - u4 * 216.5625 + u2 * 59.0625 
		- 2.1875);
	tmp3 = jr2 * 3. + jr4 * (u2 * 17.5 - 7.5) + jr6 * (u4 * 86.625 - u2 * 
		78.75 + 13.125);

	a[i__ * 3 + 1] = x[i__ * 3 + 1] * tmp1 * tmp2;
	a[i__ * 3 + 2] = x[i__ * 3 + 2] * tmp1 * tmp2;
	a[i__ * 3 + 3] = x[i__ * 3 + 3] * tmp1 * (tmp2 - tmp3);

/* Calculate barycentric accelerations on the central body */
	tmp4 = m[i__] / m[1];
	acen[1] -= tmp4 * a[i__ * 3 + 1];
	acen[2] -= tmp4 * a[i__ * 3 + 2];
	acen[3] -= tmp4 * a[i__ * 3 + 3];
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mfo_obl__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MFO_PN.FOR    (ErikSoft   3 October 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* ****** To be completed at a later date ****** */

/* Calculates post-Newtonian relativistic corrective accelerations for a set */
/* of NBOD bodies (NBIG of which are Big). */

/* This routine should not be called from the symplectic algorithm MAL_MVS */
/* or the conservative Bulirsch-Stoer algorithm MAL_BS2. */

/* N.B. All coordinates and velocities must be with respect to central body!!!! */
/* === */
/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mfo_pn__(integer *nbod, integer *nbig, doublereal *m, 
	doublereal *x, doublereal *v, doublereal *a)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;



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
    a -= 4;
    v -= 4;
    x -= 4;
    --m;

    /* Function Body */
    i__1 = *nbod;
    for (j = 1; j <= i__1; ++j) {
	a[j * 3 + 1] = 0.;
	a[j * 3 + 2] = 0.;
	a[j * 3 + 3] = 0.;
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mfo_pn__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MFO_PR.FOR    (ErikSoft   3 October 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* ****** To be completed at a later date ****** */

/* Calculates radiation pressure and Poynting-Robertson drag for a set */
/* of NBOD bodies (NBIG of which are Big). */

/* This routine should not be called from the symplectic algorithm MAL_MVS */
/* or the conservative Bulirsch-Stoer algorithm MAL_BS2. */

/* N.B. All coordinates and velocities must be with respect to central body!!!! */
/* === */
/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mfo_pr__(integer *nbod, integer *nbig, doublereal *m, 
	doublereal *x, doublereal *v, doublereal *a, doublereal *ngf)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;



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
    ngf -= 5;
    a -= 4;
    v -= 4;
    x -= 4;
    --m;

    /* Function Body */
    i__1 = *nbod;
    for (j = 1; j <= i__1; ++j) {
	a[j * 3 + 1] = 0.;
	a[j * 3 + 2] = 0.;
	a[j * 3 + 3] = 0.;
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mfo_pr__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_C2FL.FOR    (ErikSoft  1 July 1999) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Converts a CHARACTER*8 ASCII string into a REAL*8 variable. */

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

    x = mio_c2re__(c__, &c_b98, &c_b150, &c__7, (ftnlen)8);
    x = x * 2. - 1.;
    ex = *(unsigned char *)&c__[7] - 144;
    d__1 = (doublereal) ex;
    ret_val = x * pow_dd(&c_b186, &d__1);

/* ------------------------------------------------------------------------------ */

    return ret_val;
} /* mio_c2fl__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_C2RE.FOR    (ErikSoft  1 July 1999) */

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
	y = (y + (doublereal) (*(unsigned char *)&c__[j - 1] - 32)) / 224.;
    }

    ret_val = *xmin + y * (*xmax - *xmin);

/* ------------------------------------------------------------------------------ */

    return ret_val;
} /* mio_c2re__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_CE.FOR    (ErikSoft   1 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Writes details of close encounter minima to an output file, and decides how */
/* to continue the integration depending upon the close-encounter option */
/* chosen by the user. Close encounter details are stored until either 100 */
/* have been accumulated, or a data dump is done, at which point the stored */
/* encounter details are also output. */

/* For each encounter, the routine outputs the time and distance of closest */
/* approach, the identities of the objects involved, and the output */
/* variables of the objects at this time. The output variables are: */
/* expressed as */
/*  r = the radial distance */
/*  theta = polar angle */
/*  phi = azimuthal angle */
/*  fv = 1 / [1 + 2(ke/be)^2], where be and ke are the object's binding and */
/*                             kinetic energies. (Note that 0 < fv < 1). */
/*  vtheta = polar angle of velocity vector */
/*  vphi = azimuthal angle of the velocity vector */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mio_ce__(doublereal *time, doublereal *tstart, 
	doublereal *rcen, doublereal *rmax, integer *nbod, integer *nbig, 
	doublereal *m, integer *stat, char *id, integer *nclo, integer *iclo, 
	integer *jclo, integer *opt, integer *stopflag, doublereal *tclo, 
	doublereal *dclo, doublereal *ixvclo, doublereal *jxvclo, char *mem, 
	integer *lmem, char *outfile, integer *nstored, integer *ceflush, 
	ftnlen id_len, ftnlen mem_len, ftnlen outfile_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    char ch__1[8];
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    double d_lg10(doublereal *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void), f_clos(cllist *);

    /* Local variables */
    extern /* Subroutine */ int mco_x2ov__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static char c__[80*200];
    static integer k;
    static doublereal t1, fr, fv, phi, tmp0, rfac;
    static integer year;
    static doublereal vphi, theta;
    static integer month;
    static char fstop[38];
    static doublereal vtheta;
    extern /* Character */ VOID mio_fl2c__(char *, ftnlen, doublereal *);
    static char tstring[6];
    extern /* Character */ VOID mio_re2c__(char *, ftnlen, doublereal *, 
	    doublereal *, doublereal *);
    extern /* Subroutine */ int mio_jd2y__(doublereal *, integer *, integer *,
	     doublereal *);

    /* Fortran I/O blocks */
    static cilist io___685 = { 0, 22, 0, "(a1,a2,a70)", 0 };
    static cilist io___691 = { 0, 23, 0, fstop, 0 };
    static cilist io___693 = { 0, 23, 0, fstop, 0 };




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


/* Scaling factor (maximum possible range) for distances */
    /* Parameter adjustments */
    id -= 8;
    --stat;
    --m;
    jxvclo -= 7;
    ixvclo -= 7;
    --dclo;
    --tclo;
    --jclo;
    --iclo;
    --opt;
    mem -= 80;
    --lmem;
    outfile -= 80;

    /* Function Body */
    d__1 = *rmax / *rcen;
    rfac = d_lg10(&d__1);

/* Store details of each new close-encounter minimum */
    i__1 = *nclo;
    for (k = 1; k <= i__1; ++k) {
	++(*nstored);
	mio_fl2c__(ch__1, (ftnlen)8, &tclo[k]);
	s_copy(c__ + (*nstored - 1) * 80, ch__1, (ftnlen)8, (ftnlen)8);
	d__1 = (doublereal) (iclo[k] - 1);
	mio_re2c__(ch__1, (ftnlen)8, &d__1, &c_b98, &c_b99);
	s_copy(c__ + ((*nstored - 1) * 80 + 8), ch__1, (ftnlen)8, (ftnlen)8);
	d__1 = (doublereal) (jclo[k] - 1);
	mio_re2c__(ch__1, (ftnlen)8, &d__1, &c_b98, &c_b99);
	s_copy(c__ + ((*nstored - 1) * 80 + 11), ch__1, (ftnlen)8, (ftnlen)8);
	mio_fl2c__(ch__1, (ftnlen)8, &dclo[k]);
	s_copy(c__ + ((*nstored - 1) * 80 + 14), ch__1, (ftnlen)8, (ftnlen)8);

	mco_x2ov__(rcen, rmax, &m[1], &c_b98, &ixvclo[k * 6 + 1], &ixvclo[k * 
		6 + 2], &ixvclo[k * 6 + 3], &ixvclo[k * 6 + 4], &ixvclo[k * 6 
		+ 5], &ixvclo[k * 6 + 6], &fr, &theta, &phi, &fv, &vtheta, &
		vphi);
	mio_re2c__(ch__1, (ftnlen)8, &fr, &c_b98, &rfac);
	s_copy(c__ + ((*nstored - 1) * 80 + 22), ch__1, (ftnlen)8, (ftnlen)8);
	mio_re2c__(ch__1, (ftnlen)8, &theta, &c_b98, &c_b196);
	s_copy(c__ + ((*nstored - 1) * 80 + 26), ch__1, (ftnlen)8, (ftnlen)8);
	mio_re2c__(ch__1, (ftnlen)8, &phi, &c_b98, &c_b58);
	s_copy(c__ + ((*nstored - 1) * 80 + 30), ch__1, (ftnlen)8, (ftnlen)8);
	mio_re2c__(ch__1, (ftnlen)8, &fv, &c_b98, &c_b150);
	s_copy(c__ + ((*nstored - 1) * 80 + 34), ch__1, (ftnlen)8, (ftnlen)8);
	mio_re2c__(ch__1, (ftnlen)8, &vtheta, &c_b98, &c_b196);
	s_copy(c__ + ((*nstored - 1) * 80 + 38), ch__1, (ftnlen)8, (ftnlen)8);
	mio_re2c__(ch__1, (ftnlen)8, &vphi, &c_b98, &c_b58);
	s_copy(c__ + ((*nstored - 1) * 80 + 42), ch__1, (ftnlen)8, (ftnlen)8);

	mco_x2ov__(rcen, rmax, &m[1], &c_b98, &jxvclo[k * 6 + 1], &jxvclo[k * 
		6 + 2], &jxvclo[k * 6 + 3], &jxvclo[k * 6 + 4], &jxvclo[k * 6 
		+ 5], &jxvclo[k * 6 + 6], &fr, &theta, &phi, &fv, &vtheta, &
		vphi);
	mio_re2c__(ch__1, (ftnlen)8, &fr, &c_b98, &rfac);
	s_copy(c__ + ((*nstored - 1) * 80 + 46), ch__1, (ftnlen)8, (ftnlen)8);
	mio_re2c__(ch__1, (ftnlen)8, &theta, &c_b98, &c_b196);
	s_copy(c__ + ((*nstored - 1) * 80 + 50), ch__1, (ftnlen)8, (ftnlen)8);
	mio_re2c__(ch__1, (ftnlen)8, &phi, &c_b98, &c_b58);
	s_copy(c__ + ((*nstored - 1) * 80 + 54), ch__1, (ftnlen)8, (ftnlen)8);
	mio_re2c__(ch__1, (ftnlen)8, &fv, &c_b98, &c_b150);
	s_copy(c__ + ((*nstored - 1) * 80 + 58), ch__1, (ftnlen)8, (ftnlen)8);
	mio_re2c__(ch__1, (ftnlen)8, &vtheta, &c_b98, &c_b196);
	s_copy(c__ + ((*nstored - 1) * 80 + 62), ch__1, (ftnlen)12, (ftnlen)8)
		;
	mio_re2c__(ch__1, (ftnlen)8, &vphi, &c_b98, &c_b58);
	s_copy(c__ + ((*nstored - 1) * 80 + 66), ch__1, (ftnlen)12, (ftnlen)8)
		;
    }

/* If required, output the stored close encounter details */
    if (*nstored >= 100 || *ceflush == 0) {
L10:
	o__1.oerr = 1;
	o__1.ounit = 22;
	o__1.ofnmlen = 80;
	o__1.ofnm = outfile + 160;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = "append";
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L10;
	}
	i__1 = *nstored;
	for (k = 1; k <= i__1; ++k) {
	    s_wsfe(&io___685);
	    do_fio(&c__1, "\f", (ftnlen)1);
	    do_fio(&c__1, "6b", (ftnlen)2);
	    do_fio(&c__1, c__ + (k - 1) * 80, (ftnlen)70);
	    e_wsfe();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 22;
	cl__1.csta = 0;
	f_clos(&cl__1);
	*nstored = 0;
    }

/* If new encounter minima have occurred, decide whether to stop integration */
    *stopflag = 0;
    if (opt[1] == 1 && *nclo > 0) {
L20:
	o__1.oerr = 1;
	o__1.ounit = 23;
	o__1.ofnmlen = 80;
	o__1.ofnm = outfile + 240;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = "append";
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L20;
	}
/* If time style is Gregorian date then... */
	tmp0 = tclo[1];
	if (opt[3] == 1) {
	    s_copy(fstop, "(5a,/,9x,a,i10,1x,i2,1x,f4.1)", (ftnlen)38, (
		    ftnlen)29);
	    mio_jd2y__(&tmp0, &year, &month, &t1);
	    s_wsfe(&io___691);
	    do_fio(&c__1, mem + 9680, lmem[121]);
	    do_fio(&c__1, mem + 10080, lmem[126]);
	    do_fio(&c__1, id + (iclo[1] << 3), (ftnlen)8);
	    do_fio(&c__1, ",", (ftnlen)1);
	    do_fio(&c__1, id + (jclo[1] << 3), (ftnlen)8);
	    do_fio(&c__1, mem + 5680, lmem[71]);
	    do_fio(&c__1, (char *)&year, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&month, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	    e_wsfe();
/* Otherwise... */
	} else {
	    if (opt[3] == 3) {
		s_copy(tstring, mem + 160, (ftnlen)6, (ftnlen)80);
		s_copy(fstop, "(5a,/,9x,a,f14.3,a)", (ftnlen)38, (ftnlen)19);
		t1 = (tmp0 - *tstart) / 365.25;
	    } else {
		s_copy(tstring, mem + 80, (ftnlen)6, (ftnlen)80);
		s_copy(fstop, "(5a,/,9x,a,f14.1,a)", (ftnlen)38, (ftnlen)19);
		if (opt[3] == 0) {
		    t1 = tmp0;
		}
		if (opt[3] == 2) {
		    t1 = tmp0 - *tstart;
		}
	    }
	    s_wsfe(&io___693);
	    do_fio(&c__1, mem + 9680, lmem[121]);
	    do_fio(&c__1, mem + 10080, lmem[126]);
	    do_fio(&c__1, id + (iclo[1] << 3), (ftnlen)8);
	    do_fio(&c__1, ",", (ftnlen)1);
	    do_fio(&c__1, id + (jclo[1] << 3), (ftnlen)8);
	    do_fio(&c__1, mem + 5680, lmem[71]);
	    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, tstring, (ftnlen)6);
	    e_wsfe();
	}
	*stopflag = 1;
	cl__1.cerr = 0;
	cl__1.cunit = 23;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mio_ce__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_DUMP.FOR    (ErikSoft   21 February 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Writes masses, coordinates, velocities etc. of all objects, and integration */
/* parameters, to dump files. Also updates a restart file containing other */
/* variables used internally by MERCURY. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mio_dump__(doublereal *time, doublereal *tstart, 
	doublereal *tstop, doublereal *dtout, integer *algor, doublereal *h0, 
	doublereal *tol, doublereal *jcen, doublereal *rcen, doublereal *rmax,
	 doublereal *en, doublereal *am, doublereal *cefac, integer *ndump, 
	integer *nfun, integer *nbod, integer *nbig, doublereal *m, 
	doublereal *x, doublereal *v, doublereal *s, doublereal *rho, 
	doublereal *rceh, integer *stat, char *id, doublereal *ngf, 
	doublereal *epoch, integer *opt, integer *opflag, char *dumpfile, 
	char *mem, integer *lmem, ftnlen id_len, ftnlen dumpfile_len, ftnlen 
	mem_len)
{
    /* Format strings */
    static char fmt_312[] = "(1p,3(1x,e22.15),1x,i8)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;
    icilist ici__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void), s_wsle(cilist *), do_lio(integer *, 
	    integer *, char *, ftnlen), e_wsle(void), s_wsfi(icilist *), 
	    e_wsfi(void), f_clos(cllist *);

    /* Local variables */
    static char c__[150];
    static integer i__, j, k, j1, j2;
    static doublereal v0[6000]	/* was [3][2000] */, x0[6000]	/* was [3][
	    2000] */, k_2__;
    static integer idp, len;
    static doublereal rcen_2__, rcen_4__, rcen_6__, rhocgs;

    /* Fortran I/O blocks */
    static cilist io___702 = { 0, 31, 0, "(a)", 0 };
    static cilist io___705 = { 0, 31, 0, "(a)", 0 };
    static cilist io___706 = { 0, 31, 0, "(a)", 0 };
    static cilist io___707 = { 0, 31, 0, 0, 0 };
    static cilist io___708 = { 0, 31, 0, 0, 0 };
    static cilist io___709 = { 0, 31, 0, "(a)", 0 };
    static icilist io___712 = { 0, c__+8, 0, "(1p,a3,e11.5,a3,e11.5)", 29, 1 }
	    ;
    static cilist io___714 = { 0, 31, 0, "(a)", 0 };
    static cilist io___715 = { 0, 31, 0, fmt_312, 0 };
    static cilist io___717 = { 0, 31, 0, fmt_312, 0 };
    static cilist io___719 = { 0, 31, 0, fmt_312, 0 };
    static cilist io___720 = { 0, 31, 0, fmt_312, 0 };
    static cilist io___721 = { 0, 31, 0, fmt_312, 0 };
    static cilist io___722 = { 0, 33, 0, "(a)", 0 };
    static cilist io___723 = { 0, 33, 0, "(a)", 0 };
    static cilist io___724 = { 0, 33, 0, "(a)", 0 };
    static cilist io___725 = { 0, 33, 0, "(a)", 0 };
    static cilist io___726 = { 0, 33, 0, "(a)", 0 };
    static cilist io___727 = { 0, 33, 0, 0, 0 };
    static cilist io___728 = { 0, 33, 0, 0, 0 };
    static cilist io___729 = { 0, 33, 0, 0, 0 };
    static cilist io___730 = { 0, 33, 0, 0, 0 };
    static cilist io___731 = { 0, 33, 0, 0, 0 };
    static cilist io___732 = { 0, 33, 0, 0, 0 };
    static cilist io___733 = { 0, 33, 0, 0, 0 };
    static cilist io___734 = { 0, 33, 0, 0, 0 };
    static cilist io___735 = { 0, 33, 0, 0, 0 };
    static cilist io___736 = { 0, 33, 0, 0, 0 };
    static cilist io___737 = { 0, 33, 0, 0, 0 };
    static cilist io___738 = { 0, 33, 0, 0, 0 };
    static cilist io___739 = { 0, 33, 0, 0, 0 };
    static cilist io___740 = { 0, 33, 0, "(a)", 0 };
    static cilist io___741 = { 0, 33, 0, "(a)", 0 };
    static cilist io___742 = { 0, 33, 0, "(a)", 0 };
    static cilist io___743 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___744 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___745 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___746 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___747 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___748 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___749 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___750 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___751 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___752 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___753 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___754 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___755 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___756 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___757 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___758 = { 0, 33, 0, "(a)", 0 };
    static cilist io___759 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___760 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___761 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___762 = { 0, 33, 0, "(2a)", 0 };
    static cilist io___763 = { 0, 33, 0, "(a)", 0 };
    static cilist io___764 = { 0, 33, 0, "(a)", 0 };
    static cilist io___765 = { 0, 33, 0, "(a)", 0 };
    static cilist io___766 = { 0, 33, 0, 0, 0 };
    static cilist io___767 = { 0, 33, 0, 0, 0 };
    static cilist io___768 = { 0, 33, 0, 0, 0 };
    static cilist io___769 = { 0, 33, 0, 0, 0 };
    static cilist io___770 = { 0, 33, 0, 0, 0 };
    static cilist io___771 = { 0, 33, 0, 0, 0 };
    static cilist io___772 = { 0, 33, 0, 0, 0 };
    static cilist io___773 = { 0, 33, 0, 0, 0 };
    static cilist io___774 = { 0, 33, 0, 0, 0 };
    static cilist io___775 = { 0, 33, 0, 0, 0 };
    static cilist io___776 = { 0, 33, 0, 0, 0 };
    static cilist io___777 = { 0, 35, 0, "(1x,i2)", 0 };
    static cilist io___778 = { 0, 35, 0, 0, 0 };
    static cilist io___779 = { 0, 35, 0, 0, 0 };
    static cilist io___780 = { 0, 35, 0, 0, 0 };
    static cilist io___781 = { 0, 35, 0, 0, 0 };
    static cilist io___782 = { 0, 35, 0, 0, 0 };
    static cilist io___783 = { 0, 35, 0, 0, 0 };
    static cilist io___784 = { 0, 35, 0, 0, 0 };




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
    --jcen;
    --en;
    --am;
    --epoch;
    ngf -= 5;
    id -= 8;
    --stat;
    --rceh;
    --rho;
    s -= 4;
    v -= 4;
    x -= 4;
    --m;
    --opt;
    dumpfile -= 80;
    mem -= 80;
    --lmem;

    /* Function Body */
    rhocgs = 498.06095345055081;
    k_2__ = 3379.3806811609443;
    rcen_2__ = 1. / (*rcen * *rcen);
    rcen_4__ = rcen_2__ * rcen_2__;
    rcen_6__ = rcen_4__ * rcen_2__;

/* If using close-binary star, convert to user coordinates */
/*      if (algor.eq.11) call mco_h2ub (time,jcen,nbod,nbig,h0,m,x,v, */
/*     %   x0,v0) */

/* Dump to temporary files (idp=1) and real dump files (idp=2) */
    for (idp = 1; idp <= 2; ++idp) {

/* Dump data for the Big (i=1) and Small (i=2) bodies */
	for (i__ = 1; i__ <= 2; ++i__) {
	    if (idp == 1) {
		if (i__ == 1) {
		    s_copy(c__, "big.tmp     ", (ftnlen)12, (ftnlen)12);
		}
		if (i__ == 2) {
		    s_copy(c__, "small.tmp   ", (ftnlen)12, (ftnlen)12);
		}
L20:
		o__1.oerr = 1;
		o__1.ounit = 31;
		o__1.ofnmlen = 12;
		o__1.ofnm = c__;
		o__1.orl = 0;
		o__1.osta = "unknown";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		i__1 = f_open(&o__1);
		if (i__1 != 0) {
		    goto L20;
		}
	    } else {
L25:
		o__1.oerr = 1;
		o__1.ounit = 31;
		o__1.ofnmlen = 80;
		o__1.ofnm = dumpfile + i__ * 80;
		o__1.orl = 0;
		o__1.osta = "old";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		i__1 = f_open(&o__1);
		if (i__1 != 0) {
		    goto L25;
		}
	    }

/* Write header lines, data style (and epoch for Big bodies) */
	    s_wsfe(&io___702);
	    do_fio(&c__1, mem + (i__ + 151) * 80, lmem[i__ + 151]);
	    e_wsfe();
	    if (i__ == 1) {
		j1 = 2;
		j2 = *nbig;
	    } else {
		j1 = *nbig + 1;
		j2 = *nbod;
	    }
	    s_wsfe(&io___705);
	    do_fio(&c__1, mem + 12320, lmem[154]);
	    e_wsfe();
	    s_wsfe(&io___706);
	    do_fio(&c__1, mem + 12400, lmem[155]);
	    e_wsfe();
	    s_wsle(&io___707);
	    do_lio(&c__9, &c__1, mem + 12480, lmem[156]);
	    do_lio(&c__9, &c__1, "Cartesian", (ftnlen)9);
	    e_wsle();
	    if (i__ == 1) {
		s_wsle(&io___708);
		do_lio(&c__9, &c__1, mem + 12560, lmem[157]);
		do_lio(&c__5, &c__1, (char *)&(*time), (ftnlen)sizeof(
			doublereal));
		e_wsle();
	    }
	    s_wsfe(&io___709);
	    do_fio(&c__1, mem + 12400, lmem[155]);
	    e_wsfe();

/* For each body... */
	    i__1 = j2;
	    for (j = j1; j <= i__1; ++j) {
		len = 37;
		s_copy(c__, id + (j << 3), (ftnlen)8, (ftnlen)8);
		s_wsfi(&io___712);
		do_fio(&c__1, " r=", (ftnlen)3);
		do_fio(&c__1, (char *)&rceh[j], (ftnlen)sizeof(doublereal));
		do_fio(&c__1, " d=", (ftnlen)3);
		d__1 = rho[j] / rhocgs;
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		e_wsfi();
		if (m[j] > 0.) {
		    i__2 = len;
		    ici__1.icierr = 0;
		    ici__1.icirnum = 1;
		    ici__1.icirlen = len + 25 - i__2;
		    ici__1.iciunit = c__ + i__2;
		    ici__1.icifmt = "(a3,e22.15)";
		    s_wsfi(&ici__1);
		    do_fio(&c__1, " m=", (ftnlen)3);
		    d__1 = m[j] * k_2__;
		    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		    e_wsfi();
		    len += 25;
		}
		for (k = 1; k <= 3; ++k) {
		    if (ngf[k + (j << 2)] != 0.) {
			i__2 = len;
			ici__1.icierr = 0;
			ici__1.icirnum = 1;
			ici__1.icirlen = len + 16 - i__2;
			ici__1.iciunit = c__ + i__2;
			ici__1.icifmt = "(a2,i1,a1,e12.5)";
			s_wsfi(&ici__1);
			do_fio(&c__1, " a", (ftnlen)2);
			do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
			do_fio(&c__1, "=", (ftnlen)1);
			do_fio(&c__1, (char *)&ngf[k + (j << 2)], (ftnlen)
				sizeof(doublereal));
			e_wsfi();
			len += 16;
		    }
		}
		if (ngf[(j << 2) + 4] != 0.) {
		    i__2 = len;
		    ici__1.icierr = 0;
		    ici__1.icirnum = 1;
		    ici__1.icirlen = len + 15 - i__2;
		    ici__1.iciunit = c__ + i__2;
		    ici__1.icifmt = "(a3,e12.5)";
		    s_wsfi(&ici__1);
		    do_fio(&c__1, " b=", (ftnlen)3);
		    do_fio(&c__1, (char *)&ngf[(j << 2) + 4], (ftnlen)sizeof(
			    doublereal));
		    e_wsfi();
		    len += 15;
		}
		s_wsfe(&io___714);
		do_fio(&c__1, c__, len);
		e_wsfe();
		if (*algor == 11) {
		    s_wsfe(&io___715);
		    do_fio(&c__1, (char *)&x0[j * 3 - 3], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&x0[j * 3 - 2], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&x0[j * 3 - 1], (ftnlen)sizeof(
			    doublereal));
		    e_wsfe();
		    s_wsfe(&io___717);
		    do_fio(&c__1, (char *)&v0[j * 3 - 3], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&v0[j * 3 - 2], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&v0[j * 3 - 1], (ftnlen)sizeof(
			    doublereal));
		    e_wsfe();
		} else {
		    s_wsfe(&io___719);
		    do_fio(&c__1, (char *)&x[j * 3 + 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&x[j * 3 + 2], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&x[j * 3 + 3], (ftnlen)sizeof(
			    doublereal));
		    e_wsfe();
		    s_wsfe(&io___720);
		    do_fio(&c__1, (char *)&v[j * 3 + 1], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&v[j * 3 + 2], (ftnlen)sizeof(
			    doublereal));
		    do_fio(&c__1, (char *)&v[j * 3 + 3], (ftnlen)sizeof(
			    doublereal));
		    e_wsfe();
		}
		s_wsfe(&io___721);
		d__1 = s[j * 3 + 1] * k_2__;
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		d__2 = s[j * 3 + 2] * k_2__;
		do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
		d__3 = s[j * 3 + 3] * k_2__;
		do_fio(&c__1, (char *)&d__3, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = 31;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	}

/* Dump the integration parameters */
L40:
	if (idp == 1) {
	    o__1.oerr = 1;
	    o__1.ounit = 33;
	    o__1.ofnmlen = 9;
	    o__1.ofnm = "param.tmp";
	    o__1.orl = 0;
	    o__1.osta = "unknown";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L40;
	    }
	}
L45:
	if (idp == 2) {
	    o__1.oerr = 1;
	    o__1.ounit = 33;
	    o__1.ofnmlen = 80;
	    o__1.ofnm = dumpfile + 240;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L45;
	    }
	}

/* Important parameters */
	s_wsfe(&io___722);
	do_fio(&c__1, mem + 12080, lmem[151]);
	e_wsfe();
	s_wsfe(&io___723);
	do_fio(&c__1, mem + 12320, lmem[154]);
	e_wsfe();
	s_wsfe(&io___724);
	do_fio(&c__1, mem + 12400, lmem[155]);
	e_wsfe();
	s_wsfe(&io___725);
	do_fio(&c__1, mem + 12640, lmem[158]);
	e_wsfe();
	s_wsfe(&io___726);
	do_fio(&c__1, mem + 12400, lmem[155]);
	e_wsfe();
	if (*algor == 1) {
	    s_wsle(&io___727);
	    do_lio(&c__9, &c__1, mem + 12720, lmem[159]);
	    do_lio(&c__9, &c__1, "MVS", (ftnlen)3);
	    e_wsle();
	} else if (*algor == 2) {
	    s_wsle(&io___728);
	    do_lio(&c__9, &c__1, mem + 12720, lmem[159]);
	    do_lio(&c__9, &c__1, "BS", (ftnlen)2);
	    e_wsle();
	} else if (*algor == 3) {
	    s_wsle(&io___729);
	    do_lio(&c__9, &c__1, mem + 12720, lmem[159]);
	    do_lio(&c__9, &c__1, "BS2", (ftnlen)3);
	    e_wsle();
	} else if (*algor == 4) {
	    s_wsle(&io___730);
	    do_lio(&c__9, &c__1, mem + 12720, lmem[159]);
	    do_lio(&c__9, &c__1, "RADAU", (ftnlen)5);
	    e_wsle();
	} else if (*algor == 10) {
	    s_wsle(&io___731);
	    do_lio(&c__9, &c__1, mem + 12720, lmem[159]);
	    do_lio(&c__9, &c__1, "HYBRID", (ftnlen)6);
	    e_wsle();
	} else if (*algor == 11) {
	    s_wsle(&io___732);
	    do_lio(&c__9, &c__1, mem + 12720, lmem[159]);
	    do_lio(&c__9, &c__1, "CLOSE", (ftnlen)5);
	    e_wsle();
	} else if (*algor == 12) {
	    s_wsle(&io___733);
	    do_lio(&c__9, &c__1, mem + 12720, lmem[159]);
	    do_lio(&c__9, &c__1, "WIDE", (ftnlen)4);
	    e_wsle();
	} else {
	    s_wsle(&io___734);
	    do_lio(&c__9, &c__1, mem + 12720, lmem[159]);
	    do_lio(&c__9, &c__1, "0", (ftnlen)1);
	    e_wsle();
	}
	s_wsle(&io___735);
	do_lio(&c__9, &c__1, mem + 12800, lmem[160]);
	do_lio(&c__5, &c__1, (char *)&(*tstart), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___736);
	do_lio(&c__9, &c__1, mem + 12880, lmem[161]);
	do_lio(&c__5, &c__1, (char *)&(*tstop), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___737);
	do_lio(&c__9, &c__1, mem + 12960, lmem[162]);
	do_lio(&c__5, &c__1, (char *)&(*dtout), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___738);
	do_lio(&c__9, &c__1, mem + 13040, lmem[163]);
	do_lio(&c__5, &c__1, (char *)&(*h0), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___739);
	do_lio(&c__9, &c__1, mem + 13120, lmem[164]);
	do_lio(&c__5, &c__1, (char *)&(*tol), (ftnlen)sizeof(doublereal));
	e_wsle();

/* Integration options */
	s_wsfe(&io___740);
	do_fio(&c__1, mem + 12400, lmem[155]);
	e_wsfe();
	s_wsfe(&io___741);
	do_fio(&c__1, mem + 13200, lmem[165]);
	e_wsfe();
	s_wsfe(&io___742);
	do_fio(&c__1, mem + 12400, lmem[155]);
	e_wsfe();
	if (opt[1] == 0) {
	    s_wsfe(&io___743);
	    do_fio(&c__1, mem + 13280, lmem[166]);
	    do_fio(&c__1, mem + 400, lmem[5]);
	    e_wsfe();
	} else {
	    s_wsfe(&io___744);
	    do_fio(&c__1, mem + 13280, lmem[166]);
	    do_fio(&c__1, mem + 480, lmem[6]);
	    e_wsfe();
	}
	if (opt[2] == 0) {
	    s_wsfe(&io___745);
	    do_fio(&c__1, mem + 13360, lmem[167]);
	    do_fio(&c__1, mem + 400, lmem[5]);
	    e_wsfe();
	    s_wsfe(&io___746);
	    do_fio(&c__1, mem + 13440, lmem[168]);
	    do_fio(&c__1, mem + 400, lmem[5]);
	    e_wsfe();
	} else if (opt[2] == 2) {
	    s_wsfe(&io___747);
	    do_fio(&c__1, mem + 13360, lmem[167]);
	    do_fio(&c__1, mem + 480, lmem[6]);
	    e_wsfe();
	    s_wsfe(&io___748);
	    do_fio(&c__1, mem + 13440, lmem[168]);
	    do_fio(&c__1, mem + 480, lmem[6]);
	    e_wsfe();
	} else {
	    s_wsfe(&io___749);
	    do_fio(&c__1, mem + 13360, lmem[167]);
	    do_fio(&c__1, mem + 480, lmem[6]);
	    e_wsfe();
	    s_wsfe(&io___750);
	    do_fio(&c__1, mem + 13440, lmem[168]);
	    do_fio(&c__1, mem + 400, lmem[5]);
	    e_wsfe();
	}
	if (opt[3] == 0 || opt[3] == 2) {
	    s_wsfe(&io___751);
	    do_fio(&c__1, mem + 13520, lmem[169]);
	    do_fio(&c__1, mem + 80, lmem[1]);
	    e_wsfe();
	} else {
	    s_wsfe(&io___752);
	    do_fio(&c__1, mem + 13520, lmem[169]);
	    do_fio(&c__1, mem + 160, lmem[2]);
	    e_wsfe();
	}
	if (opt[3] == 2 || opt[3] == 3) {
	    s_wsfe(&io___753);
	    do_fio(&c__1, mem + 13600, lmem[170]);
	    do_fio(&c__1, mem + 480, lmem[6]);
	    e_wsfe();
	} else {
	    s_wsfe(&io___754);
	    do_fio(&c__1, mem + 13600, lmem[170]);
	    do_fio(&c__1, mem + 400, lmem[5]);
	    e_wsfe();
	}
	if (opt[4] == 1) {
	    s_wsfe(&io___755);
	    do_fio(&c__1, mem + 13680, lmem[171]);
	    do_fio(&c__1, mem + 560, lmem[7]);
	    e_wsfe();
	} else if (opt[4] == 3) {
	    s_wsfe(&io___756);
	    do_fio(&c__1, mem + 13680, lmem[171]);
	    do_fio(&c__1, mem + 720, lmem[9]);
	    e_wsfe();
	} else {
	    s_wsfe(&io___757);
	    do_fio(&c__1, mem + 13680, lmem[171]);
	    do_fio(&c__1, mem + 640, lmem[8]);
	    e_wsfe();
	}
	s_wsfe(&io___758);
	do_fio(&c__1, mem + 13760, lmem[172]);
	e_wsfe();
	if (opt[7] == 1) {
	    s_wsfe(&io___759);
	    do_fio(&c__1, mem + 13840, lmem[173]);
	    do_fio(&c__1, mem + 480, lmem[6]);
	    e_wsfe();
	} else {
	    s_wsfe(&io___760);
	    do_fio(&c__1, mem + 13840, lmem[173]);
	    do_fio(&c__1, mem + 400, lmem[5]);
	    e_wsfe();
	}
	if (opt[8] == 1) {
	    s_wsfe(&io___761);
	    do_fio(&c__1, mem + 13920, lmem[174]);
	    do_fio(&c__1, mem + 480, lmem[6]);
	    e_wsfe();
	} else {
	    s_wsfe(&io___762);
	    do_fio(&c__1, mem + 13920, lmem[174]);
	    do_fio(&c__1, mem + 400, lmem[5]);
	    e_wsfe();
	}

/* Infrequently-changed parameters */
	s_wsfe(&io___763);
	do_fio(&c__1, mem + 12400, lmem[155]);
	e_wsfe();
	s_wsfe(&io___764);
	do_fio(&c__1, mem + 14000, lmem[175]);
	e_wsfe();
	s_wsfe(&io___765);
	do_fio(&c__1, mem + 12400, lmem[155]);
	e_wsfe();
	s_wsle(&io___766);
	do_lio(&c__9, &c__1, mem + 14080, lmem[176]);
	do_lio(&c__5, &c__1, (char *)&(*rmax), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___767);
	do_lio(&c__9, &c__1, mem + 14160, lmem[177]);
	do_lio(&c__5, &c__1, (char *)&(*rcen), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___768);
	do_lio(&c__9, &c__1, mem + 14240, lmem[178]);
	d__1 = m[1] * k_2__;
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___769);
	do_lio(&c__9, &c__1, mem + 14320, lmem[179]);
	d__1 = jcen[1] * rcen_2__;
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___770);
	do_lio(&c__9, &c__1, mem + 14400, lmem[180]);
	d__1 = jcen[2] * rcen_4__;
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___771);
	do_lio(&c__9, &c__1, mem + 14480, lmem[181]);
	d__1 = jcen[3] * rcen_6__;
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___772);
	do_lio(&c__9, &c__1, mem + 14560, lmem[182]);
	e_wsle();
	s_wsle(&io___773);
	do_lio(&c__9, &c__1, mem + 14640, lmem[183]);
	e_wsle();
	s_wsle(&io___774);
	do_lio(&c__9, &c__1, mem + 14720, lmem[184]);
	do_lio(&c__5, &c__1, (char *)&(*cefac), (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___775);
	do_lio(&c__9, &c__1, mem + 14800, lmem[185]);
	do_lio(&c__3, &c__1, (char *)&(*ndump), (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___776);
	do_lio(&c__9, &c__1, mem + 14880, lmem[186]);
	do_lio(&c__3, &c__1, (char *)&(*nfun), (ftnlen)sizeof(integer));
	e_wsle();
	cl__1.cerr = 0;
	cl__1.cunit = 33;
	cl__1.csta = 0;
	f_clos(&cl__1);

/* Create new version of the restart file */
L60:
	if (idp == 1) {
	    o__1.oerr = 1;
	    o__1.ounit = 35;
	    o__1.ofnmlen = 11;
	    o__1.ofnm = "restart.tmp";
	    o__1.orl = 0;
	    o__1.osta = "unknown";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L60;
	    }
	}
L65:
	if (idp == 2) {
	    o__1.oerr = 1;
	    o__1.ounit = 35;
	    o__1.ofnmlen = 80;
	    o__1.ofnm = dumpfile + 320;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L65;
	    }
	}
	s_wsfe(&io___777);
	do_fio(&c__1, (char *)&(*opflag), (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsle(&io___778);
	d__1 = en[1] * k_2__;
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___779);
	d__1 = am[1] * k_2__;
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___780);
	d__1 = en[3] * k_2__;
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___781);
	d__1 = am[3] * k_2__;
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___782);
	d__1 = s[4] * k_2__;
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___783);
	d__1 = s[5] * k_2__;
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
	s_wsle(&io___784);
	d__1 = s[6] * k_2__;
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsle();
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/* ------------------------------------------------------------------------------ */

/* L311: */
/* L313: */
/* L314: */
    return 0;
} /* mio_dump__ */


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
    static cilist io___785 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___786 = { 0, 0, 0, "(/,3a,/,2a)", 0 };




/* Input/Output */

/* ------------------------------------------------------------------------------ */

    s_wsfe(&io___785);
    do_fio(&c__1, " ERROR: Programme terminated. See information", (ftnlen)45)
	    ;
    do_fio(&c__1, " file for details.", (ftnlen)18);
    e_wsfe();

    io___786.ciunit = *unit;
    s_wsfe(&io___786);
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

/*      MIO_FL2C.FOR    (ErikSoft  1 July 1998) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Converts a (floating point) REAL*8 variable X, into a CHARACTER*8 ASCII */
/* string, using the new format compression: */

/* X is first converted to base 224, and then each base 224 digit is converted */
/* to an ASCII character, such that 0 -> character 32, 1 -> character 33... */
/* and 223 -> character 255. */
/* The first 7 characters in the string are used to store the mantissa, and the */
/* eighth character is used for the exponent. */

/* ASCII characters 0 - 31 (CTRL characters) are not used, because they */
/* cause problems when using some operating systems. */

/* N.B. X must lie in the range -1.e112 < X < 1.e112 */
/* === */

/* ------------------------------------------------------------------------------ */

/* Character */ VOID mio_fl2c__(char *ret_val, ftnlen ret_val_len, doublereal 
	*x)
{
    /* System generated locals */
    integer i__1;
    char ch__2[8];

    /* Builtin functions */
    double d_lg10(doublereal *), pow_di(doublereal *, integer *), d_sign(
	    doublereal *, doublereal *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static doublereal y, ax;
    static integer ex;
    extern /* Character */ VOID mio_re2c__(char *, ftnlen, doublereal *, 
	    doublereal *, doublereal *);



/* Input/Output */

/* Local */

/* ------------------------------------------------------------------------------ */

    if (*x == 0.) {
	y = .5;
    } else {
	ax = abs(*x);
	ex = (integer) d_lg10(&ax);
	if (ax >= 1.) {
	    ++ex;
	}
	i__1 = -ex;
	y = ax * pow_di(&c_b186, &i__1);
	if (y == 1.) {
	    y *= .1;
	    ++ex;
	}
	y = d_sign(&y, x) * .5 + .5;
    }

    mio_re2c__(ch__2, (ftnlen)8, &y, &c_b98, &c_b150);
    s_copy(ret_val, ch__2, (ftnlen)8, (ftnlen)8);
    ex += 112;
    if (ex > 223) {
	ex = 223;
    }
    if (ex < 0) {
	ex = 0;
    }
    *(unsigned char *)&ret_val[7] = (char) (ex + 32);

/* ------------------------------------------------------------------------------ */

    return ;
} /* mio_fl2c__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_IN.FOR    (ErikSoft   4 May 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Reads names, masses, coordinates and velocities of all the bodies, */
/* and integration parameters for the MERCURY integrator package. */
/* If DUMPFILE(4) exists, the routine assumes this is a continuation of */
/* an old integration, and reads all the data from the dump files instead */
/* of the input files. */

/* N.B. All coordinates are with respect to the central body!! */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mio_in__(doublereal *time, doublereal *tstart, 
	doublereal *tstop, doublereal *dtout, integer *algor, doublereal *h0, 
	doublereal *tol, doublereal *rmax, doublereal *rcen, doublereal *jcen,
	 doublereal *en, doublereal *am, doublereal *cefac, integer *ndump, 
	integer *nfun, integer *nbod, integer *nbig, doublereal *m, 
	doublereal *x, doublereal *v, doublereal *s, doublereal *rho, 
	doublereal *rceh, integer *stat, char *id, doublereal *epoch, 
	doublereal *ngf, integer *opt, integer *opflag, integer *ngflag, char 
	*outfile, char *dumpfile, integer *lmem, char *mem, ftnlen id_len, 
	ftnlen outfile_len, ftnlen dumpfile_len, ftnlen mem_len)
{
    /* Initialized data */

    static char alg[3*60] = "MVS" "Mvs" "mvs" "mvs" "mvs" "BS " "Bs " "bs " 
	    "Bul" "bul" "BS2" "Bs2" "bs2" "Bu2" "bu2" "RAD" "Rad" "rad" "RA " 
	    "ra " "xxx" "xxx" "xxx" "xxx" "xxx" "xxx" "xxx" "xxx" "xxx" "xxx" 
	    "xxx" "xxx" "xxx" "xxx" "xxx" "xxx" "xxx" "xxx" "xxx" "xxx" "TES" 
	    "Tes" "tes" "Tst" "tst" "HYB" "Hyb" "hyb" "HY " "hy " "CLO" "Clo" 
	    "clo" "CB " "cb " "WID" "Wid" "wid" "WB " "wb ";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    icilist ici__1;
    olist o__1;
    cllist cl__1;
    alist al__1;
    inlist ioin__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_inqu(inlist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer f_open(olist *), s_rsfe(cilist *), e_rsfe(void), f_clos(cllist *),
	     s_cmp(char *, char *, ftnlen, ftnlen), s_rsli(icilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsli(void);
    double pow_dd(doublereal *, doublereal *);
    integer f_back(alist *), s_rsle(cilist *), e_rsle(void);
    double sqrt(doublereal), d_mod(doublereal *, doublereal *);
    integer s_wsfi(icilist *), e_wsfi(void);

    /* Local variables */
    static char filename[80];
    static integer informat;
    static doublereal a, e, i__;
    static integer j, k;
    static doublereal l, n, p, q;
    static char c1[1], c3[3];
    static doublereal t1;
    static char c80[80];
    static integer lim[20]	/* was [2][10] */;
    static doublereal tmp2, tmp3, tmp4, tmp5, tmp6;
    static integer year;
    static doublereal temp;
    static integer nsub, itmp, jtmp;
    static logical test, flag1, flag2;
    static integer month;
    static char infile__[80*3];
    static integer lineno;
    static doublereal rhocgs;
    extern /* Subroutine */ int mxx_en__(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static char string[150];
    static logical oldflag;
    extern /* Subroutine */ int mio_err__(integer *, char *, integer *, char *
	    , integer *, char *, integer *, char *, integer *, ftnlen, ftnlen,
	     ftnlen, ftnlen), mio_spl__(integer *, char *, integer *, integer 
	    *, ftnlen), mco_el2x__(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), mco_x2el__(doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), mio_jd2y__(doublereal *
	    , integer *, integer *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___796 = { 0, 6, 0, "(/,2a)", 0 };
    static cilist io___797 = { 0, 16, 1, "(i3,1x,i2,1x,a80)", 0 };
    static cilist io___798 = { 0, 15, 0, "(a150)", 0 };
    static cilist io___803 = { 0, 15, 0, "(a150)", 0 };
    static cilist io___804 = { 0, 15, 0, "(a150)", 0 };
    static cilist io___807 = { 0, 13, 0, "(a150)", 0 };
    static icilist io___809 = { 1, c80, 0, 0, 80, 1 };
    static icilist io___810 = { 1, c80, 0, 0, 80, 1 };
    static icilist io___811 = { 1, c80, 0, 0, 80, 1 };
    static icilist io___812 = { 1, c80, 0, 0, 80, 1 };
    static icilist io___813 = { 1, c80, 0, 0, 80, 1 };
    static icilist io___815 = { 1, c80, 0, 0, 80, 1 };
    static icilist io___816 = { 1, c80, 0, 0, 80, 1 };
    static icilist io___817 = { 1, c80, 0, 0, 80, 1 };
    static icilist io___818 = { 1, c80, 0, 0, 80, 1 };
    static icilist io___819 = { 1, c80, 0, 0, 80, 1 };
    static icilist io___820 = { 1, c80, 0, 0, 80, 1 };
    static icilist io___821 = { 1, c80, 0, 0, 80, 1 };
    static icilist io___822 = { 1, c80, 0, 0, 80, 1 };
    static icilist io___823 = { 1, c80, 0, 0, 80, 1 };
    static cilist io___825 = { 0, 13, 0, "(/,2a)", 0 };
    static cilist io___826 = { 0, 11, 0, "(a150)", 0 };
    static cilist io___829 = { 0, 11, 0, "(a150)", 0 };
    static cilist io___830 = { 0, 11, 1, "(a)", 0 };
    static cilist io___831 = { 0, 23, 0, "(/,3a)", 0 };
    static cilist io___833 = { 0, 11, 1, "(a150)", 0 };
    static cilist io___834 = { 1, 11, 0, 0, 0 };
    static cilist io___835 = { 1, 11, 0, 0, 0 };
    static cilist io___846 = { 0, 23, 0, "(/,a,i10,i2,f8.5,/)", 0 };
    static cilist io___847 = { 0, 23, 0, "(/,a,f18.7,a,/)", 0 };
    static cilist io___848 = { 0, 23, 0, "(/,a,f18.5,a,/)", 0 };
    static cilist io___849 = { 0, 35, 0, 0, 0 };
    static cilist io___850 = { 0, 35, 0, 0, 0 };
    static cilist io___851 = { 0, 35, 0, 0, 0 };
    static cilist io___852 = { 0, 23, 0, "(/,a)", 0 };
    static cilist io___853 = { 0, 23, 0, "(a)", 0 };
    static cilist io___854 = { 0, 23, 0, "(/,2a)", 0 };
    static cilist io___855 = { 0, 23, 0, "(/,a,1p,e19.13,a)", 0 };
    static cilist io___856 = { 0, 23, 0, "(/,a,f19.7,a)", 0 };
    static cilist io___857 = { 0, 23, 0, "(a,1p,e19.13)", 0 };
    static cilist io___858 = { 0, 23, 0, "(a,f19.7)", 0 };
    static cilist io___859 = { 0, 23, 0, "(a,f15.3)", 0 };
    static cilist io___860 = { 0, 23, 0, "(2a)", 0 };
    static cilist io___861 = { 0, 23, 0, "(2a)", 0 };
    static cilist io___862 = { 0, 23, 0, "(2a)", 0 };
    static cilist io___863 = { 0, 23, 0, "(/,a,f10.3,a)", 0 };
    static cilist io___864 = { 0, 23, 0, "(a,1p1e10.4)", 0 };
    static cilist io___865 = { 0, 23, 0, "(a,1p1e10.4,a)", 0 };
    static cilist io___866 = { 0, 23, 0, "(a,1p1e11.4)", 0 };
    static cilist io___867 = { 0, 23, 0, "(a,1p1e11.4)", 0 };
    static cilist io___868 = { 0, 23, 0, "(a,1p1e11.4)", 0 };
    static cilist io___869 = { 0, 23, 0, "(a,1p1e10.4,a)", 0 };
    static cilist io___870 = { 0, 23, 0, "(a,1p1e10.4,a)", 0 };
    static cilist io___872 = { 0, 23, 0, "(/,2a)", 0 };
    static cilist io___873 = { 0, 23, 0, "(2a)", 0 };
    static cilist io___874 = { 0, 23, 0, "(2a)", 0 };
    static cilist io___875 = { 0, 23, 0, "(2a)", 0 };
    static cilist io___876 = { 0, 23, 0, "(/,a,i4)", 0 };
    static cilist io___877 = { 0, 23, 0, "(a,i4)", 0 };
    static cilist io___878 = { 0, 23, 0, "(//,a)", 0 };
    static cilist io___879 = { 0, 23, 0, "(a)", 0 };
    static cilist io___880 = { 0, 23, 0, "(/,a,1p1e12.5,a)", 0 };
    static cilist io___881 = { 0, 23, 0, "(a,1p1e12.5,a)", 0 };
    static cilist io___882 = { 0, 23, 0, "(/,2a)", 0 };
    static cilist io___883 = { 0, 23, 0, "(/,2a)", 0 };
    static cilist io___884 = { 0, 23, 0, "(/,2a)", 0 };
    static cilist io___892 = { 0, 23, 0, "(/,2a,/a)", 0 };
    static icilist io___893 = { 0, c3, 0, "(i3)", 3, 1 };




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
/*      real*8 v0(3,NMAX),x0(3,NMAX) */

/* ------------------------------------------------------------------------------ */

    /* Parameter adjustments */
    mem -= 80;
    --lmem;
    dumpfile -= 80;
    outfile -= 80;
    --opt;
    ngf -= 5;
    --epoch;
    id -= 8;
    --stat;
    --rceh;
    --rho;
    s -= 4;
    v -= 4;
    x -= 4;
    --m;
    --am;
    --en;
    --jcen;

    /* Function Body */

    rhocgs = 498.06095345055081;
    for (j = 1; j <= 80; ++j) {
	*(unsigned char *)&filename[j - 1] = ' ';
    }
    for (j = 1; j <= 3; ++j) {
	s_copy(infile__ + (j - 1) * 80, filename, (ftnlen)80, (ftnlen)80);
	s_copy(outfile + j * 80, filename, (ftnlen)80, (ftnlen)80);
	s_copy(dumpfile + j * 80, filename, (ftnlen)80, (ftnlen)80);
    }
    s_copy(dumpfile + 320, filename, (ftnlen)80, (ftnlen)80);

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
	s_wsfe(&io___796);
	do_fio(&c__1, " ERROR: This file is needed to start", (ftnlen)36);
	do_fio(&c__1, " the integration:  message.in", (ftnlen)29);
	e_wsfe();
	s_stop("", (ftnlen)0);
    }
    o__1.oerr = 0;
    o__1.ounit = 16;
    o__1.ofnmlen = 10;
    o__1.ofnm = "message.in";
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
L10:
    i__1 = s_rsfe(&io___797);
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = do_fio(&c__1, (char *)&lmem[j], (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L20;
    }
    i__1 = do_fio(&c__1, mem + j * 80, (ftnlen)80);
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
    cl__1.cunit = 16;
    cl__1.csta = 0;
    f_clos(&cl__1);

/* Read in filenames and check for duplicate filenames */
    ioin__1.inerr = 0;
    ioin__1.infilen = 8;
    ioin__1.infile = "files.in";
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
	mio_err__(&c__6, mem + 6480, &lmem[81], mem + 7040, &lmem[88], " ", &
		c__1, "files.in", &c__8, (ftnlen)80, (ftnlen)80, (ftnlen)1, (
		ftnlen)8);
    }
    o__1.oerr = 0;
    o__1.ounit = 15;
    o__1.ofnmlen = 8;
    o__1.ofnm = "files.in";
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);

/* Input files */
    for (j = 1; j <= 3; ++j) {
	s_rsfe(&io___798);
	do_fio(&c__1, string, (ftnlen)150);
	e_rsfe();
	mio_spl__(&c__150, string, &nsub, lim, (ftnlen)150);
	i__1 = lim[0] - 1;
	s_copy(infile__ + (j - 1) * 80, string + i__1, lim[1] - lim[0] + 1, 
		lim[1] - i__1);
	i__1 = j - 1;
	for (k = 1; k <= i__1; ++k) {
	    if (s_cmp(infile__ + (j - 1) * 80, infile__ + (k - 1) * 80, (
		    ftnlen)80, (ftnlen)80) == 0) {
		mio_err__(&c__6, mem + 6480, &lmem[81], mem + 7120, &lmem[89],
			 infile__ + (j - 1) * 80, &c__80, mem + 6880, &lmem[
			86], (ftnlen)80, (ftnlen)80, (ftnlen)80, (ftnlen)80);
	    }
	}
    }

/* Output files */
    for (j = 1; j <= 3; ++j) {
	s_rsfe(&io___803);
	do_fio(&c__1, string, (ftnlen)150);
	e_rsfe();
	mio_spl__(&c__150, string, &nsub, lim, (ftnlen)150);
	i__1 = lim[0] - 1;
	s_copy(outfile + j * 80, string + i__1, lim[1] - lim[0] + 1, lim[1] - 
		i__1);
	i__1 = j - 1;
	for (k = 1; k <= i__1; ++k) {
	    if (s_cmp(outfile + j * 80, outfile + k * 80, (ftnlen)80, (ftnlen)
		    80) == 0) {
		mio_err__(&c__6, mem + 6480, &lmem[81], mem + 7120, &lmem[89],
			 outfile + j * 80, &c__80, mem + 6880, &lmem[86], (
			ftnlen)80, (ftnlen)80, (ftnlen)80, (ftnlen)80);
	    }
	}
	for (k = 1; k <= 3; ++k) {
	    if (s_cmp(outfile + j * 80, infile__ + (k - 1) * 80, (ftnlen)80, (
		    ftnlen)80) == 0) {
		mio_err__(&c__6, mem + 6480, &lmem[81], mem + 7120, &lmem[89],
			 outfile + j * 80, &c__80, mem + 6880, &lmem[86], (
			ftnlen)80, (ftnlen)80, (ftnlen)80, (ftnlen)80);
	    }
	}
    }

/* Dump files */
    for (j = 1; j <= 4; ++j) {
	s_rsfe(&io___804);
	do_fio(&c__1, string, (ftnlen)150);
	e_rsfe();
	mio_spl__(&c__150, string, &nsub, lim, (ftnlen)150);
	i__1 = lim[0] - 1;
	s_copy(dumpfile + j * 80, string + i__1, lim[1] - lim[0] + 1, lim[1] 
		- i__1);
	i__1 = j - 1;
	for (k = 1; k <= i__1; ++k) {
	    if (s_cmp(dumpfile + j * 80, dumpfile + k * 80, (ftnlen)80, (
		    ftnlen)80) == 0) {
		mio_err__(&c__6, mem + 6480, &lmem[81], mem + 7120, &lmem[89],
			 dumpfile + j * 80, &c__80, mem + 6880, &lmem[86], (
			ftnlen)80, (ftnlen)80, (ftnlen)80, (ftnlen)80);
	    }
	}
	for (k = 1; k <= 3; ++k) {
	    if (s_cmp(dumpfile + j * 80, infile__ + (k - 1) * 80, (ftnlen)80, 
		    (ftnlen)80) == 0) {
		mio_err__(&c__6, mem + 6480, &lmem[81], mem + 7120, &lmem[89],
			 dumpfile + j * 80, &c__80, mem + 6880, &lmem[86], (
			ftnlen)80, (ftnlen)80, (ftnlen)80, (ftnlen)80);
	    }
	}
	for (k = 1; k <= 3; ++k) {
	    if (s_cmp(dumpfile + j * 80, outfile + k * 80, (ftnlen)80, (
		    ftnlen)80) == 0) {
		mio_err__(&c__6, mem + 6480, &lmem[81], mem + 7120, &lmem[89],
			 dumpfile + j * 80, &c__80, mem + 6880, &lmem[86], (
			ftnlen)80, (ftnlen)80, (ftnlen)80, (ftnlen)80);
	    }
	}
    }
    cl__1.cerr = 0;
    cl__1.cunit = 15;
    cl__1.csta = 0;
    f_clos(&cl__1);

/* Find out if this is an old integration (i.e. does the restart file exist) */
    ioin__1.inerr = 0;
    ioin__1.infilen = 80;
    ioin__1.infile = dumpfile + 320;
    ioin__1.inex = &oldflag;
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

/* Check if information file exists, and append a continuation message */
    if (oldflag) {
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = outfile + 240;
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
	    mio_err__(&c__6, mem + 6480, &lmem[81], mem + 7040, &lmem[88], 
		    " ", &c__1, outfile + 240, &c__80, (ftnlen)80, (ftnlen)80,
		     (ftnlen)1, (ftnlen)80);
	}
L320:
	o__1.oerr = 1;
	o__1.ounit = 23;
	o__1.ofnmlen = 80;
	o__1.ofnm = outfile + 240;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = "append";
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L320;
	}
    } else {

/* If new integration, check information file doesn't exist, and then create it */
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = outfile + 240;
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
	    mio_err__(&c__6, mem + 6480, &lmem[81], mem + 6960, &lmem[87], 
		    " ", &c__1, outfile + 240, &c__80, (ftnlen)80, (ftnlen)80,
		     (ftnlen)1, (ftnlen)80);
	}
L410:
	o__1.oerr = 1;
	o__1.ounit = 23;
	o__1.ofnmlen = 80;
	o__1.ofnm = outfile + 240;
	o__1.orl = 0;
	o__1.osta = "new";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L410;
	}
    }

/* ------------------------------------------------------------------------------ */

/*  READ  IN  INTEGRATION  PARAMETERS */

/* Check if the file containing integration parameters exists, and open it */
    s_copy(filename, infile__ + 160, (ftnlen)80, (ftnlen)80);
    if (oldflag) {
	s_copy(filename, dumpfile + 240, (ftnlen)80, (ftnlen)80);
    }
    ioin__1.inerr = 0;
    ioin__1.infilen = 80;
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
    if (! test) {
	mio_err__(&c__23, mem + 6480, &lmem[81], mem + 7040, &lmem[88], " ", &
		c__1, filename, &c__80, (ftnlen)80, (ftnlen)80, (ftnlen)1, (
		ftnlen)80);
    }
L30:
    o__1.oerr = 1;
    o__1.ounit = 13;
    o__1.ofnmlen = 80;
    o__1.ofnm = filename;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L30;
    }

/* Read integration parameters */
    lineno = 0;
    for (j = 1; j <= 26; ++j) {
L40:
	++lineno;
	s_rsfe(&io___807);
	do_fio(&c__1, string, (ftnlen)150);
	e_rsfe();
	if (*(unsigned char *)string == ')') {
	    goto L40;
	}
	mio_spl__(&c__150, string, &nsub, lim, (ftnlen)150);
	s_copy(c80, "   ", (ftnlen)3, (ftnlen)3);
	i__1 = lim[(nsub << 1) - 2] - 1;
	s_copy(c80, string + i__1, (ftnlen)80, lim[(nsub << 1) - 1] - i__1);
	if (j == 1) {
	    *algor = 0;
	    for (k = 1; k <= 60; ++k) {
		if (s_cmp(c80, alg + (k - 1) * 3, (ftnlen)3, (ftnlen)3) == 0) 
			{
		    *algor = (k + 4) / 5;
		}
	    }
	    if (*algor == 0) {
		i__1 = lim[(nsub << 1) - 2] - 1;
		i__2 = lim[(nsub << 1) - 1] - lim[(nsub << 1) - 2] + 1;
		mio_err__(&c__23, mem + 6480, &lmem[81], mem + 7840, &lmem[98]
			, c80 + i__1, &i__2, mem + 6800, &lmem[85], (ftnlen)
			80, (ftnlen)80, lim[(nsub << 1) - 1] - i__1, (ftnlen)
			80);
	    }
	}
	if (j == 2) {
	    i__1 = s_rsli(&io___809);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*tstart), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
	if (j == 3) {
	    i__1 = s_rsli(&io___810);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*tstop), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
	if (j == 4) {
	    i__1 = s_rsli(&io___811);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*dtout), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
	if (j == 5) {
	    i__1 = s_rsli(&io___812);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*h0), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
	if (j == 6) {
	    i__1 = s_rsli(&io___813);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*tol), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
	*(unsigned char *)c1 = *(unsigned char *)c80;
	if (j == 7 && (*(unsigned char *)c1 == 'y' || *(unsigned char *)c1 == 
		'Y')) {
	    opt[1] = 1;
	}
	if (j == 8 && (*(unsigned char *)c1 == 'n' || *(unsigned char *)c1 == 
		'N')) {
	    opt[2] = 0;
	}
	if (j == 9 && (*(unsigned char *)c1 == 'y' || *(unsigned char *)c1 == 
		'Y')) {
	    opt[2] = 2;
	}
	if (j == 10 && (*(unsigned char *)c1 == 'd' || *(unsigned char *)c1 ==
		 'D')) {
	    opt[3] = 0;
	}
	if (j == 11 && (*(unsigned char *)c1 == 'y' || *(unsigned char *)c1 ==
		 'Y')) {
	    opt[3] += 2;
	}
	if (j == 12) {
	    if (*(unsigned char *)c1 == 'l' || *(unsigned char *)c1 == 'L') {
		opt[4] = 1;
	    } else if (j == 12 && (*(unsigned char *)c1 == 'm' || *(unsigned 
		    char *)c1 == 'M')) {
		opt[4] = 2;
	    } else if (j == 12 && (*(unsigned char *)c1 == 'h' || *(unsigned 
		    char *)c1 == 'H')) {
		opt[4] = 3;
	    } else {
		goto L661;
	    }
	}
	if (j == 15 && (*(unsigned char *)c1 == 'y' || *(unsigned char *)c1 ==
		 'Y')) {
	    opt[8] = 1;
	}
	if (j == 16) {
	    i__1 = s_rsli(&io___815);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*rmax), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
	if (j == 17) {
	    i__1 = s_rsli(&io___816);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*rcen), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
	if (j == 18) {
	    i__1 = s_rsli(&io___817);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&m[1], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
	if (j == 19) {
	    i__1 = s_rsli(&io___818);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&jcen[1], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
	if (j == 20) {
	    i__1 = s_rsli(&io___819);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&jcen[2], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
	if (j == 21) {
	    i__1 = s_rsli(&io___820);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&jcen[3], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
	if (j == 24) {
	    i__1 = s_rsli(&io___821);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*cefac), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
	if (j == 25) {
	    i__1 = s_rsli(&io___822);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&(*ndump), (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
	if (j == 26) {
	    i__1 = s_rsli(&io___823);
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = do_lio(&c__3, &c__1, (char *)&(*nfun), (ftnlen)sizeof(
		    integer));
	    if (i__1 != 0) {
		goto L661;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L661;
	    }
	}
    }
    *h0 = abs(*h0);
    *tol = abs(*tol);
    *rmax = abs(*rmax);
    *rcen = abs(*rcen);
    *cefac = abs(*cefac);
    cl__1.cerr = 0;
    cl__1.cunit = 13;
    cl__1.csta = 0;
    f_clos(&cl__1);

/* Change quantities for central object to suitable units */
    m[1] = abs(m[1]) * 2.959122082855911e-4;
    jcen[1] = jcen[1] * *rcen * *rcen;
    jcen[2] = jcen[2] * *rcen * *rcen * *rcen * *rcen;
    jcen[3] = jcen[3] * *rcen * *rcen * *rcen * *rcen * *rcen * *rcen;
    s[4] = 0.;
    s[5] = 0.;
    s[6] = 0.;

/* Make sure that RCEN isn't too small, since it is used to scale the output */
/* (Minimum value corresponds to a central body with density 100g/cm^3). */
    temp = pow_dd(&m[1], &c_b666) * .0011235;
    if (*rcen < temp) {
	*rcen = temp;
	s_wsfe(&io___825);
	do_fio(&c__1, mem + 9680, lmem[121]);
	do_fio(&c__1, mem + 10480, lmem[131]);
	e_wsfe();
    }

/* ------------------------------------------------------------------------------ */

/*  READ  IN  DATA  FOR  BIG  AND  SMALL  BODIES */

    *nbod = 1;
    for (j = 1; j <= 2; ++j) {
	if (j == 2) {
	    *nbig = *nbod;
	}

/* Check if the file containing data for Big bodies exists, and open it */
	s_copy(filename, infile__ + (j - 1) * 80, (ftnlen)80, (ftnlen)80);
	if (oldflag) {
	    s_copy(filename, dumpfile + j * 80, (ftnlen)80, (ftnlen)80);
	}
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
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
	if (! test) {
	    mio_err__(&c__23, mem + 6480, &lmem[81], mem + 7040, &lmem[88], 
		    " ", &c__1, filename, &c__80, (ftnlen)80, (ftnlen)80, (
		    ftnlen)1, (ftnlen)80);
	}
L110:
	o__1.oerr = 1;
	o__1.ounit = 11;
	o__1.ofnmlen = 80;
	o__1.ofnm = filename;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L110;
	}

/* Read data style */
L120:
	s_rsfe(&io___826);
	do_fio(&c__1, string, (ftnlen)150);
	e_rsfe();
	if (*(unsigned char *)string == ')') {
	    goto L120;
	}
	mio_spl__(&c__150, string, &nsub, lim, (ftnlen)150);
	i__1 = lim[(nsub << 1) - 2] - 1;
	s_copy(c3, string + i__1, (ftnlen)3, lim[(nsub << 1) - 2] + 2 - i__1);
	if (s_cmp(c3, "Car", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "car", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "CAR", (ftnlen)3, (
		ftnlen)3) == 0) {
	    informat = 1;
	} else if (s_cmp(c3, "Ast", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, 
		"ast", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "AST", (ftnlen)
		3, (ftnlen)3) == 0) {
	    informat = 2;
	} else if (s_cmp(c3, "Com", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, 
		"com", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "COM", (ftnlen)
		3, (ftnlen)3) == 0) {
	    informat = 3;
	} else {
	    mio_err__(&c__23, mem + 6480, &lmem[81], mem + 7280, &lmem[91], 
		    " ", &c__1, mem + (j + 82) * 80, &lmem[j + 82], (ftnlen)
		    80, (ftnlen)80, (ftnlen)1, (ftnlen)80);
	}

/* Read epoch of Big bodies */
	if (j == 1) {
L125:
	    s_rsfe(&io___829);
	    do_fio(&c__1, string, (ftnlen)150);
	    e_rsfe();
	    if (*(unsigned char *)string == ')') {
		goto L125;
	    }
	    mio_spl__(&c__150, string, &nsub, lim, (ftnlen)150);
	    i__2 = lim[(nsub << 1) - 2] - 1;
	    ici__1.icierr = 1;
	    ici__1.iciend = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = lim[(nsub << 1) - 1] - i__2;
	    ici__1.iciunit = string + i__2;
	    ici__1.icifmt = 0;
	    i__1 = s_rsli(&ici__1);
	    if (i__1 != 0) {
		goto L667;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&(*time), (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L667;
	    }
	    i__1 = e_rsli();
	    if (i__1 != 0) {
		goto L667;
	    }
	}

/* Read information for each object */
L130:
	i__1 = s_rsfe(&io___830);
	if (i__1 != 0) {
	    goto L140;
	}
	i__1 = do_fio(&c__1, string, (ftnlen)150);
	if (i__1 != 0) {
	    goto L140;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L140;
	}
	if (*(unsigned char *)string == ')') {
	    goto L130;
	}
	mio_spl__(&c__150, string, &nsub, lim, (ftnlen)150);
	if (lim[0] == -1) {
	    goto L140;
	}

/* Determine the name of the object */
	++(*nbod);
	if (*nbod > 2000) {
	    mio_err__(&c__23, mem + 6480, &lmem[81], mem + 7200, &lmem[90], 
		    " ", &c__1, mem + 6560, &lmem[82], (ftnlen)80, (ftnlen)80,
		     (ftnlen)1, (ftnlen)80);
	}

	if (lim[1] - lim[0] > 7) {
	    s_wsfe(&io___831);
	    do_fio(&c__1, mem + 9680, lmem[121]);
	    do_fio(&c__1, mem + 9760, lmem[122]);
	    i__1 = lim[0] - 1;
	    do_fio(&c__1, string + i__1, lim[1] - i__1);
	    e_wsfe();
	}
	i__1 = lim[0] - 1;
/* Computing MIN */
	i__2 = lim[0] + 7;
	s_copy(id + (*nbod << 3), string + i__1, (ftnlen)8, min(i__2,lim[1]) 
		- i__1);
/* Check if another object has the same name */
	i__1 = *nbod - 1;
	for (k = 1; k <= i__1; ++k) {
	    if (s_cmp(id + (k << 3), id + (*nbod << 3), (ftnlen)8, (ftnlen)8) 
		    == 0) {
		mio_err__(&c__23, mem + 6480, &lmem[81], mem + 8240, &lmem[
			103], id + (*nbod << 3), &c__8, " ", &c__1, (ftnlen)
			80, (ftnlen)80, (ftnlen)8, (ftnlen)1);
	    }
	}

/* Default values of mass, close-encounter limit, density etc. */
	m[*nbod] = 0.;
	rceh[*nbod] = 1.;
	rho[*nbod] = rhocgs;
	epoch[*nbod] = *time;
	for (k = 1; k <= 4; ++k) {
	    ngf[k + (*nbod << 2)] = 0.;
	}

/* Read values of mass, close-encounter limit, density etc. */
	i__1 = nsub;
	for (k = 3; k <= i__1; k += 2) {
	    i__2 = lim[(k - 1 << 1) - 2] - 1;
	    s_copy(c80, string + i__2, (ftnlen)80, lim[(k - 1 << 1) - 1] - 
		    i__2);
	    i__3 = lim[(k << 1) - 2] - 1;
	    ici__1.icierr = 1;
	    ici__1.iciend = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = lim[(k << 1) - 1] - i__3;
	    ici__1.iciunit = string + i__3;
	    ici__1.icifmt = 0;
	    i__2 = s_rsli(&ici__1);
	    if (i__2 != 0) {
		goto L666;
	    }
	    i__2 = do_lio(&c__5, &c__1, (char *)&temp, (ftnlen)sizeof(
		    doublereal));
	    if (i__2 != 0) {
		goto L666;
	    }
	    i__2 = e_rsli();
	    if (i__2 != 0) {
		goto L666;
	    }
	    if (*(unsigned char *)c80 == 'm' || *(unsigned char *)c80 == 'M') 
		    {
		m[*nbod] = temp * 2.959122082855911e-4;
	    } else if (*(unsigned char *)c80 == 'r' || *(unsigned char *)c80 
		    == 'R') {
		rceh[*nbod] = temp;
	    } else if (*(unsigned char *)c80 == 'd' || *(unsigned char *)c80 
		    == 'D') {
		rho[*nbod] = temp * rhocgs;
	    } else if (m[*nbod] < 0. || rceh[*nbod] < 0. || rho[*nbod] < 0.) {
		mio_err__(&c__23, mem + 6480, &lmem[81], mem + 7760, &lmem[97]
			, id + (*nbod << 3), &c__8, mem + (j + 82) * 80, &
			lmem[j + 82], (ftnlen)80, (ftnlen)80, (ftnlen)8, (
			ftnlen)80);
	    } else if (s_cmp(c80, "ep", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c80, "EP", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c80, "Ep", 
		    (ftnlen)2, (ftnlen)2) == 0) {
		epoch[*nbod] = temp;
	    } else if (s_cmp(c80, "a1", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c80, "A1", (ftnlen)2, (ftnlen)2) == 0) {
		ngf[(*nbod << 2) + 1] = temp;
	    } else if (s_cmp(c80, "a2", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c80, "A2", (ftnlen)2, (ftnlen)2) == 0) {
		ngf[(*nbod << 2) + 2] = temp;
	    } else if (s_cmp(c80, "a3", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c80, "A3", (ftnlen)2, (ftnlen)2) == 0) {
		ngf[(*nbod << 2) + 3] = temp;
	    } else if (*(unsigned char *)c80 == 'b' || *(unsigned char *)c80 
		    == 'B') {
		ngf[(*nbod << 2) + 4] = temp;
	    } else {
		goto L666;
	    }
	}

/* If required, read Cartesian coordinates, velocities and spins of the bodies */
	jtmp = 100;
L135:
	i__1 = s_rsfe(&io___833);
	if (i__1 != 0) {
	    goto L666;
	}
	i__1 = do_fio(&c__1, string, (ftnlen)150);
	if (i__1 != 0) {
	    goto L666;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L666;
	}
	if (*(unsigned char *)string == ')') {
	    goto L135;
	}
	al__1.aerr = 0;
	al__1.aunit = 11;
	f_back(&al__1);
	if (informat == 1) {
	    i__1 = s_rsle(&io___834);
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&x[*nbod * 3 + 1], (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&x[*nbod * 3 + 2], (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&x[*nbod * 3 + 3], (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&v[*nbod * 3 + 1], (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&v[*nbod * 3 + 2], (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&v[*nbod * 3 + 3], (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&s[*nbod * 3 + 1], (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&s[*nbod * 3 + 2], (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&s[*nbod * 3 + 3], (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L666;
	    }
	} else {
	    i__1 = s_rsle(&io___835);
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&a, (ftnlen)sizeof(doublereal)
		    );
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&e, (ftnlen)sizeof(doublereal)
		    );
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&i__, (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&p, (ftnlen)sizeof(doublereal)
		    );
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&n, (ftnlen)sizeof(doublereal)
		    );
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&l, (ftnlen)sizeof(doublereal)
		    );
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&s[*nbod * 3 + 1], (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&s[*nbod * 3 + 2], (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = do_lio(&c__5, &c__1, (char *)&s[*nbod * 3 + 3], (ftnlen)
		    sizeof(doublereal));
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__1 = e_rsle();
	    if (i__1 != 0) {
		goto L666;
	    }
	    i__ *= .017453292519943295;
	    p = (p + n) * .017453292519943295;
	    n *= .017453292519943295;
	    temp = m[*nbod] + m[1];

/* Alternatively, read Cometary or asteroidal elements */
	    if (informat == 3) {
		q = a;
		a = q / (1. - e);
		d__2 = sqrt(temp / (d__1 = a * a * a, abs(d__1))) * (epoch[*
			nbod] - l);
		l = d_mod(&d__2, &c_b58);
	    } else {
		q = a * (1. - e);
		l *= .017453292519943295;
	    }
	    if (*algor == 11 && *nbod != 2) {
		temp += m[2];
	    }
	    mco_el2x__(&temp, &q, &e, &i__, &p, &n, &l, &x[*nbod * 3 + 1], &x[
		    *nbod * 3 + 2], &x[*nbod * 3 + 3], &v[*nbod * 3 + 1], &v[*
		    nbod * 3 + 2], &v[*nbod * 3 + 3]);
	}

	s[*nbod * 3 + 1] *= 2.959122082855911e-4;
	s[*nbod * 3 + 2] *= 2.959122082855911e-4;
	s[*nbod * 3 + 3] *= 2.959122082855911e-4;

	goto L130;
L140:
	cl__1.cerr = 0;
	cl__1.cunit = 11;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

/* Set non-gravitational-forces flag, NGFLAG */
    *ngflag = 0;
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	if (ngf[(j << 2) + 1] != 0. || ngf[(j << 2) + 2] != 0. || ngf[(j << 2)
		 + 3] != 0.) {
	    if (*ngflag == 0) {
		*ngflag = 1;
	    }
	    if (*ngflag == 2) {
		*ngflag = 3;
	    }
	} else if (ngf[(j << 2) + 4] != 0.) {
	    if (*ngflag == 0) {
		*ngflag = 2;
	    }
	    if (*ngflag == 1) {
		*ngflag = 3;
	    }
	}
    }

/* ------------------------------------------------------------------------------ */

/*  IF  CONTINUING  AN  OLD  INTEGRATION */

    if (oldflag) {
	if (opt[3] == 1) {
	    mio_jd2y__(time, &year, &month, &t1);
	    s_wsfe(&io___846);
	    do_fio(&c__1, mem + 4960, lmem[62]);
	    do_fio(&c__1, (char *)&year, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&month, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else if (opt[3] == 3) {
	    t1 = (*time - *tstart) / 365.25;
	    s_wsfe(&io___847);
	    do_fio(&c__1, mem + 4960, lmem[62]);
	    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, mem + 160, lmem[2]);
	    e_wsfe();
	} else {
	    if (opt[3] == 0) {
		t1 = *time;
	    }
	    if (opt[3] == 2) {
		t1 = *time - *tstart;
	    }
	    s_wsfe(&io___848);
	    do_fio(&c__1, mem + 4960, lmem[62]);
	    do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, mem + 80, lmem[1]);
	    e_wsfe();
	}

/* Read in energy and angular momentum variables, and convert to internal units */
L330:
	o__1.oerr = 1;
	o__1.ounit = 35;
	o__1.ofnmlen = 80;
	o__1.ofnm = dumpfile + 320;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L330;
	}
	s_rsle(&io___849);
	do_lio(&c__3, &c__1, (char *)&(*opflag), (ftnlen)sizeof(integer));
	e_rsle();
	s_rsle(&io___850);
	do_lio(&c__5, &c__1, (char *)&en[1], (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&am[1], (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&en[3], (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&am[3], (ftnlen)sizeof(doublereal));
	e_rsle();
	en[1] *= 2.959122082855911e-4;
	en[3] *= 2.959122082855911e-4;
	am[1] *= 2.959122082855911e-4;
	am[3] *= 2.959122082855911e-4;
	s_rsle(&io___851);
	do_lio(&c__5, &c__1, (char *)&s[4], (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&s[5], (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&s[6], (ftnlen)sizeof(doublereal));
	e_rsle();
	s[4] *= 2.959122082855911e-4;
	s[5] *= 2.959122082855911e-4;
	s[6] *= 2.959122082855911e-4;
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
	if (*opflag == 0) {
	    *opflag = 1;
	}

/* ------------------------------------------------------------------------------ */

/*  IF  STARTING  A  NEW  INTEGRATION */

    } else {
	*opflag = -2;

/* Write integration parameters to information file */
	s_wsfe(&io___852);
	do_fio(&c__1, mem + 880, lmem[11]);
	e_wsfe();
	s_wsfe(&io___853);
	do_fio(&c__1, mem + 960, lmem[12]);
	e_wsfe();
	j = *algor + 13;
	s_wsfe(&io___854);
	do_fio(&c__1, mem + 1040, lmem[13]);
	do_fio(&c__1, mem + j * 80, lmem[j]);
	e_wsfe();
	if (*tstart >= 1e11 || *tstart <= -1e10) {
	    s_wsfe(&io___855);
	    do_fio(&c__1, mem + 2080, lmem[26]);
	    do_fio(&c__1, (char *)&(*tstart), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, mem + 80, lmem[1]);
	    e_wsfe();
	} else {
	    s_wsfe(&io___856);
	    do_fio(&c__1, mem + 2080, lmem[26]);
	    do_fio(&c__1, (char *)&(*tstart), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, mem + 80, lmem[1]);
	    e_wsfe();
	}
	if (*tstop >= 1e11 || *tstop <= -1e10) {
	    s_wsfe(&io___857);
	    do_fio(&c__1, mem + 2160, lmem[27]);
	    do_fio(&c__1, (char *)&(*tstop), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    s_wsfe(&io___858);
	    do_fio(&c__1, mem + 2160, lmem[27]);
	    do_fio(&c__1, (char *)&(*tstop), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
	s_wsfe(&io___859);
	do_fio(&c__1, mem + 2240, lmem[28]);
	do_fio(&c__1, (char *)&(*dtout), (ftnlen)sizeof(doublereal));
	e_wsfe();
	if (opt[4] == 1) {
	    s_wsfe(&io___860);
	    do_fio(&c__1, mem + 3200, lmem[40]);
	    do_fio(&c__1, mem + 560, lmem[7]);
	    e_wsfe();
	}
	if (opt[4] == 2) {
	    s_wsfe(&io___861);
	    do_fio(&c__1, mem + 3200, lmem[40]);
	    do_fio(&c__1, mem + 640, lmem[8]);
	    e_wsfe();
	}
	if (opt[4] == 3) {
	    s_wsfe(&io___862);
	    do_fio(&c__1, mem + 3200, lmem[40]);
	    do_fio(&c__1, mem + 720, lmem[9]);
	    e_wsfe();
	}

	s_wsfe(&io___863);
	do_fio(&c__1, mem + 2400, lmem[30]);
	do_fio(&c__1, (char *)&(*h0), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, mem + 80, lmem[1]);
	e_wsfe();
	s_wsfe(&io___864);
	do_fio(&c__1, mem + 2480, lmem[31]);
	do_fio(&c__1, (char *)&(*tol), (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsfe(&io___865);
	do_fio(&c__1, mem + 2560, lmem[32]);
	d__1 = m[1] / 2.959122082855911e-4;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, mem + 240, lmem[3]);
	e_wsfe();
	s_wsfe(&io___866);
	do_fio(&c__1, mem + 2640, lmem[33]);
/* Computing 2nd power */
	d__2 = *rcen;
	d__1 = jcen[1] / (d__2 * d__2);
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsfe(&io___867);
	do_fio(&c__1, mem + 2720, lmem[34]);
/* Computing 4th power */
	d__2 = *rcen, d__2 *= d__2;
	d__1 = jcen[2] / (d__2 * d__2);
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsfe(&io___868);
	do_fio(&c__1, mem + 2800, lmem[35]);
/* Computing 6th power */
	d__2 = *rcen, d__2 *= d__2;
	d__1 = jcen[3] / (d__2 * (d__2 * d__2));
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsfe(&io___869);
	do_fio(&c__1, mem + 2880, lmem[36]);
	do_fio(&c__1, (char *)&(*rmax), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, mem + 320, lmem[4]);
	e_wsfe();
	s_wsfe(&io___870);
	do_fio(&c__1, mem + 2960, lmem[37]);
	do_fio(&c__1, (char *)&(*rcen), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, mem + 320, lmem[4]);
	e_wsfe();

	itmp = 5;
	if (opt[2] == 1 || opt[2] == 2) {
	    itmp = 6;
	}
	s_wsfe(&io___872);
	do_fio(&c__1, mem + 3280, lmem[41]);
	do_fio(&c__1, mem + itmp * 80, lmem[itmp]);
	e_wsfe();
	itmp = 5;
	if (opt[2] == 2) {
	    itmp = 6;
	}
	s_wsfe(&io___873);
	do_fio(&c__1, mem + 3360, lmem[42]);
	do_fio(&c__1, mem + itmp * 80, lmem[itmp]);
	e_wsfe();
	itmp = 5;
	if (opt[7] == 1) {
	    itmp = 6;
	}
	s_wsfe(&io___874);
	do_fio(&c__1, mem + 3600, lmem[45]);
	do_fio(&c__1, mem + itmp * 80, lmem[itmp]);
	e_wsfe();
	itmp = 5;
	if (opt[8] == 1) {
	    itmp = 6;
	}
	s_wsfe(&io___875);
	do_fio(&c__1, mem + 3680, lmem[46]);
	do_fio(&c__1, mem + itmp * 80, lmem[itmp]);
	e_wsfe();

/* Check that element and close-encounter files don't exist, and create them */
	for (j = 1; j <= 2; ++j) {
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 80;
	    ioin__1.infile = outfile + j * 80;
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
		mio_err__(&c__23, mem + 6480, &lmem[81], mem + 6960, &lmem[87]
			, " ", &c__1, outfile + j * 80, &c__80, (ftnlen)80, (
			ftnlen)80, (ftnlen)1, (ftnlen)80);
	    }
L430:
	    o__1.oerr = 1;
	    o__1.ounit = j + 20;
	    o__1.ofnmlen = 80;
	    o__1.ofnm = outfile + j * 80;
	    o__1.orl = 0;
	    o__1.osta = "new";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L430;
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = j + 20;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	}

/* Check that dump files don't exist, and then create them */
	for (j = 1; j <= 4; ++j) {
	    ioin__1.inerr = 0;
	    ioin__1.infilen = 80;
	    ioin__1.infile = dumpfile + j * 80;
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
		mio_err__(&c__23, mem + 6480, &lmem[81], mem + 6960, &lmem[87]
			, " ", &c__1, dumpfile + j * 80, &c__80, (ftnlen)80, (
			ftnlen)80, (ftnlen)1, (ftnlen)80);
	    }
L450:
	    o__1.oerr = 1;
	    o__1.ounit = j + 30;
	    o__1.ofnmlen = 80;
	    o__1.ofnm = dumpfile + j * 80;
	    o__1.orl = 0;
	    o__1.osta = "new";
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L450;
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = j + 30;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	}

/* Write number of Big bodies and Small bodies to information file */
	s_wsfe(&io___876);
	do_fio(&c__1, mem + 3040, lmem[38]);
	i__1 = *nbig - 1;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___877);
	do_fio(&c__1, mem + 3120, lmem[39]);
	i__1 = *nbod - *nbig;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	e_wsfe();

/* Calculate initial energy and angular momentum and write to information file */
	s[4] = 0.;
	s[5] = 0.;
	s[6] = 0.;
	mxx_en__(&jcen[1], nbod, nbig, &m[1], &x[4], &v[4], &s[4], &en[1], &
		am[1]);
	s_wsfe(&io___878);
	do_fio(&c__1, mem + 4080, lmem[51]);
	e_wsfe();
	s_wsfe(&io___879);
	do_fio(&c__1, mem + 4160, lmem[52]);
	e_wsfe();
	s_wsfe(&io___880);
	do_fio(&c__1, mem + 4240, lmem[53]);
	d__1 = en[1] / 2.959122082855911e-4;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, mem + 5760, lmem[72]);
	e_wsfe();
	s_wsfe(&io___881);
	do_fio(&c__1, mem + 4320, lmem[54]);
	d__1 = am[1] / 2.959122082855911e-4;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, mem + 5840, lmem[73]);
	e_wsfe();

/* Initialize lost energy and angular momentum */
	en[3] = 0.;
	am[3] = 0.;

/* Write warning messages if necessary */
	if (*tstop < *tstart) {
	    s_wsfe(&io___882);
	    do_fio(&c__1, mem + 9680, lmem[121]);
	    do_fio(&c__1, mem + 9840, lmem[123]);
	    e_wsfe();
	}
	if (*nbig <= 0) {
	    s_wsfe(&io___883);
	    do_fio(&c__1, mem + 9680, lmem[121]);
	    do_fio(&c__1, mem + 9920, lmem[124]);
	    e_wsfe();
	}
	if (*nbig == *nbod) {
	    s_wsfe(&io___884);
	    do_fio(&c__1, mem + 9680, lmem[121]);
	    do_fio(&c__1, mem + 10000, lmem[125]);
	    e_wsfe();
	}
    }

/* ------------------------------------------------------------------------------ */

/*  CHECK  FOR  ATTEMPTS  TO  DO  INCOMPATIBLE  THINGS */

/* If using close-binary algorithm, set radius of central body to be no less */
/* than the periastron of binary star. */
    if (*algor == 11) {
	temp = m[1] + m[2];
	mco_x2el__(&temp, &x[7], &x[8], &x[9], &v[7], &v[8], &v[9], &a, &tmp2,
		 &tmp3, &tmp4, &tmp5, &tmp6);
	*rcen = max(*rcen,a);
    }

/* Check if non-grav forces are being used with an incompatible algorithm */
    if (*ngflag != 0 && (*algor == 3 || *algor == 11 || *algor == 12)) {
	mio_err__(&c__23, mem + 6480, &lmem[81], mem + 7360, &lmem[92], " ", &
		c__1, mem + 6800, &lmem[85], (ftnlen)80, (ftnlen)80, (ftnlen)
		1, (ftnlen)80);
    }

/* Check if user-defined force routine is being used with wrong algorithm */
    if (opt[8] == 1 && (*algor == 11 || *algor == 12)) {
	mio_err__(&c__23, mem + 6480, &lmem[81], mem + 7440, &lmem[93], " ", &
		c__1, mem + 6800, &lmem[85], (ftnlen)80, (ftnlen)80, (ftnlen)
		1, (ftnlen)80);
    }

/* Check whether MVS is being used to integrate massive Small bodies, */
/* or whether massive Small bodies have different epochs than Big bodies. */
    flag1 = FALSE_;
    flag2 = FALSE_;
    i__1 = *nbod;
    for (j = *nbig + 1; j <= i__1; ++j) {
	if (m[j] != 0.) {
	    if (*algor == 1) {
		mio_err__(&c__23, mem + 6480, &lmem[81], mem + 7520, &lmem[94]
			, " ", &c__1, mem + 6800, &lmem[85], (ftnlen)80, (
			ftnlen)80, (ftnlen)1, (ftnlen)80);
	    }
	    flag1 = TRUE_;
	}
	if (epoch[j] != *time) {
	    flag2 = TRUE_;
	}
    }
    if (flag1 && flag2) {
	mio_err__(&c__23, mem + 6480, &lmem[81], mem + 7600, &lmem[95], " ", &
		c__1, mem + 6720, &lmem[84], (ftnlen)80, (ftnlen)80, (ftnlen)
		1, (ftnlen)80);
    }

/* Check if central oblateness is being used with close-binary algorithm */
    if (*algor == 11 && (jcen[1] != 0. || jcen[2] != 0. || jcen[3] != 0.)) {
	mio_err__(&c__23, mem + 6480, &lmem[81], mem + 8160, &lmem[102], 
		" ", &c__1, mem + 6800, &lmem[85], (ftnlen)80, (ftnlen)80, (
		ftnlen)1, (ftnlen)80);
    }

/* Check whether RCEN > RMAX or RMAX/RCEN is very large */
    if (*rcen > *rmax) {
	mio_err__(&c__23, mem + 6480, &lmem[81], mem + 8400, &lmem[105], 
		" ", &c__1, mem + 6800, &lmem[85], (ftnlen)80, (ftnlen)80, (
		ftnlen)1, (ftnlen)80);
    }
    if (*rmax / *rcen >= 1e12) {
	s_wsfe(&io___892);
	do_fio(&c__1, mem + 9680, lmem[121]);
	do_fio(&c__1, mem + 8480, lmem[106]);
	do_fio(&c__1, mem + 6800, lmem[85]);
	e_wsfe();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 23;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;

/* Error reading from the input file containing integration parameters */
L661:
    s_wsfi(&io___893);
    do_fio(&c__1, (char *)&lineno, (ftnlen)sizeof(integer));
    e_wsfi();
    mio_err__(&c__23, mem + 6480, &lmem[81], mem + 7920, &lmem[99], c3, &c__3,
	     mem + 6800, &lmem[85], (ftnlen)80, (ftnlen)80, (ftnlen)3, (
	    ftnlen)80);

/* Error reading from the input file for Big or Small bodies */
L666:
    mio_err__(&c__23, mem + 6480, &lmem[81], mem + 8000, &lmem[100], id + (*
	    nbod << 3), &c__8, mem + (j + 82) * 80, &lmem[j + 82], (ftnlen)80,
	     (ftnlen)80, (ftnlen)8, (ftnlen)80);

/* Error reading epoch of Big bodies */
L667:
    mio_err__(&c__23, mem + 6480, &lmem[81], mem + 8080, &lmem[101], " ", &
	    c__1, mem + 6640, &lmem[83], (ftnlen)80, (ftnlen)80, (ftnlen)1, (
	    ftnlen)80);

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mio_in__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_JD2Y.FOR    (ErikSoft  7 July 1999) */

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

/* Subroutine */ int mio_jd2y__(doublereal *jd0, integer *year, integer *
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

/* Algorithm for negative Julian day numbers (Duffett-Smith doesn't work) */
    x = *jd0 - 2232101.5f;
    f = x - d_int(&x);
    if (f < 0.) {
	f += 1.;
    }
    d__1 = d_mod(&x, &c_b957) + 1461.;
    y = d_int(&d__1);
    d__1 = d_mod(&y, &c_b958);
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
} /* mio_jd2y__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_LOG.FOR    (ErikSoft   25 February 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Writes a progress report to the log file (or the screen if you are running */
/* Mercury interactively). */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mio_log__(doublereal *time, doublereal *tstart, 
	doublereal *en, doublereal *am, integer *opt, char *mem, integer *
	lmem, ftnlen mem_len)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal t1, tmp0, tmp1;
    static char flog[38];
    static integer year, month;
    static char tstring[6];
    extern /* Subroutine */ int mio_jd2y__(doublereal *, integer *, integer *,
	     doublereal *);

    /* Fortran I/O blocks */
    static cilist io___914 = { 0, 6, 0, flog, 0 };
    static cilist io___915 = { 0, 6, 0, flog, 0 };




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
    --opt;
    --am;
    --en;

    /* Function Body */
    if (opt[3] == 0 || opt[3] == 2) {
	s_copy(tstring, mem + 80, (ftnlen)6, (ftnlen)80);
	s_copy(flog, "(1x,a,f14.1,a,2(a,1p1e12.5))", (ftnlen)38, (ftnlen)28);
    } else if (opt[3] == 1) {
	s_copy(flog, "(1x,a,i10,1x,i2,1x,f4.1,2(a,1p1e12.5))", (ftnlen)38, (
		ftnlen)38);
    } else {
	s_copy(tstring, mem + 160, (ftnlen)6, (ftnlen)80);
	s_copy(flog, "(1x,a,f14.3,a,2(a,1p1e12.5))", (ftnlen)38, (ftnlen)28);
    }

    tmp0 = 0.;
    tmp1 = 0.;
    if (en[1] != 0.) {
	tmp0 = (en[2] + en[3] - en[1]) / abs(en[1]);
    }
    if (am[1] != 0.) {
	tmp1 = (am[2] + am[3] - am[1]) / abs(am[1]);
    }

    if (opt[3] == 1) {
	mio_jd2y__(time, &year, &month, &t1);
	s_wsfe(&io___914);
	do_fio(&c__1, mem + 5120, lmem[64]);
	do_fio(&c__1, (char *)&year, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&month, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, mem + 5200, lmem[65]);
	do_fio(&c__1, (char *)&tmp0, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, mem + 5280, lmem[66]);
	do_fio(&c__1, (char *)&tmp1, (ftnlen)sizeof(doublereal));
	e_wsfe();
    } else {
	if (opt[3] == 0) {
	    t1 = *time;
	}
	if (opt[3] == 2) {
	    t1 = *time - *tstart;
	}
	if (opt[3] == 3) {
	    t1 = (*time - *tstart) / 365.25;
	}
	s_wsfe(&io___915);
	do_fio(&c__1, mem + 5040, lmem[63]);
	do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, tstring, (ftnlen)6);
	do_fio(&c__1, mem + 5200, lmem[65]);
	do_fio(&c__1, (char *)&tmp0, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, mem + 5280, lmem[66]);
	do_fio(&c__1, (char *)&tmp1, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mio_log__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_OUT.FOR    (ErikSoft   13 February 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Writes output variables for each object to an output file. Each variable */
/* is scaled between the minimum and maximum possible values and then */
/* written in a compressed format using ASCII characters. */
/* The output variables are: */
/*  r = the radial distance */
/*  theta = polar angle */
/*  phi = azimuthal angle */
/*  fv = 1 / [1 + 2(ke/be)^2], where be and ke are the object's binding and */
/*                             kinetic energies. (Note that 0 < fv < 1). */
/*  vtheta = polar angle of velocity vector */
/*  vphi = azimuthal angle of the velocity vector */

/* If this is the first output (OPFLAG = -1), or the first output since the */
/* number of the objects or their masses have changed (OPFLAG = 1), then */
/* the names, masses and spin components of all the objects are also output. */

/* N.B. Each object's distance must lie between RCEN < R < RMAX */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mio_out__(doublereal *time, doublereal *jcen, doublereal 
	*rcen, doublereal *rmax, integer *nbod, integer *nbig, doublereal *m, 
	doublereal *xh, doublereal *vh, doublereal *s, doublereal *rho, 
	integer *stat, char *id, integer *opt, integer *opflag, integer *
	algor, char *outfile, ftnlen id_len, ftnlen outfile_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    char ch__1[8];
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    double d_lg10(doublereal *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , f_open(olist *), s_wsfe(cilist *), e_wsfe(void), f_clos(cllist *
	    );

    /* Local variables */
    extern /* Subroutine */ int mco_x2ov__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static char c__[80*2000];
    static integer k;
    static doublereal fr, fv, k_2__;
    static integer len;
    static doublereal phi, rfac, vphi;
    static char fout[5];
    static integer nchar;
    static doublereal theta, rcen_2__;
    static char header[80];
    static doublereal rhocgs, vtheta;
    extern /* Character */ VOID mio_fl2c__(char *, ftnlen, doublereal *), 
	    mio_re2c__(char *, ftnlen, doublereal *, doublereal *, doublereal 
	    *);

    /* Fortran I/O blocks */
    static icilist io___923 = { 0, fout+2, 0, "(i1)", 1, 1 };
    static icilist io___924 = { 0, fout+2, 0, "(i2)", 2, 1 };
    static cilist io___928 = { 0, 21, 0, "(a1,a2,i2,a62,i1)", 0 };
    static cilist io___929 = { 0, 21, 0, "(a51)", 0 };
    static cilist io___936 = { 0, 21, 0, "(a1,a2,a14)", 0 };
    static cilist io___937 = { 0, 21, 0, fout, 0 };




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
    --jcen;
    id -= 8;
    --stat;
    --rho;
    s -= 4;
    vh -= 4;
    xh -= 4;
    --m;
    --opt;

    /* Function Body */
    rhocgs = 498.06095345055081;
    k_2__ = 3379.3806811609443;
    rcen_2__ = 1. / (*rcen * *rcen);

/* Scaling factor (maximum possible range) for distances */
    d__1 = *rmax / *rcen;
    rfac = d_lg10(&d__1);

/* Create the format list, FOUT, used when outputting the orbital elements */
    if (opt[4] == 1) {
	nchar = 2;
    }
    if (opt[4] == 2) {
	nchar = 4;
    }
    if (opt[4] == 3) {
	nchar = 7;
    }
    len = nchar * 6 + 3;
    s_copy(fout, "(a  )", (ftnlen)5, (ftnlen)5);
    if (len < 10) {
	s_wsfi(&io___923);
	do_fio(&c__1, (char *)&len, (ftnlen)sizeof(integer));
	e_wsfi();
    }
    if (len >= 10) {
	s_wsfi(&io___924);
	do_fio(&c__1, (char *)&len, (ftnlen)sizeof(integer));
	e_wsfi();
    }

/* Open the orbital-elements output file */
L10:
    o__1.oerr = 1;
    o__1.ounit = 21;
    o__1.ofnmlen = 80;
    o__1.ofnm = outfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = "append";
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L10;
    }

/* ------------------------------------------------------------------------------ */

/*  SPECIAL  OUTPUT  PROCEDURE */

/* If this is a new integration or a complete output is required (e.g. because */
/* the number of objects has changed), then output object details & parameters. */
    if (*opflag == -1 || *opflag == 1) {

/* Compose a header line with time, number of objects and relevant parameters */
	mio_fl2c__(ch__1, (ftnlen)8, time);
	s_copy(header, ch__1, (ftnlen)8, (ftnlen)8);
	d__1 = (doublereal) (*nbig - 1);
	mio_re2c__(ch__1, (ftnlen)8, &d__1, &c_b98, &c_b99);
	s_copy(header + 8, ch__1, (ftnlen)8, (ftnlen)8);
	d__1 = (doublereal) (*nbod - *nbig);
	mio_re2c__(ch__1, (ftnlen)8, &d__1, &c_b98, &c_b99);
	s_copy(header + 11, ch__1, (ftnlen)8, (ftnlen)8);
	d__1 = m[1] * k_2__;
	mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	s_copy(header + 14, ch__1, (ftnlen)8, (ftnlen)8);
	d__1 = jcen[1] * rcen_2__;
	mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	s_copy(header + 22, ch__1, (ftnlen)8, (ftnlen)8);
	d__1 = jcen[2] * rcen_2__ * rcen_2__;
	mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	s_copy(header + 30, ch__1, (ftnlen)8, (ftnlen)8);
	d__1 = jcen[3] * rcen_2__ * rcen_2__ * rcen_2__;
	mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	s_copy(header + 38, ch__1, (ftnlen)8, (ftnlen)8);
	mio_fl2c__(ch__1, (ftnlen)8, rcen);
	s_copy(header + 46, ch__1, (ftnlen)8, (ftnlen)8);
	mio_fl2c__(ch__1, (ftnlen)8, rmax);
	s_copy(header + 54, ch__1, (ftnlen)8, (ftnlen)8);

/* For each object, compress its index number, name, mass, spin components */
/* and density (some of these need to be converted to normal units). */
	i__1 = *nbod;
	for (k = 2; k <= i__1; ++k) {
	    d__1 = (doublereal) (k - 1);
	    mio_re2c__(ch__1, (ftnlen)8, &d__1, &c_b98, &c_b99);
	    s_copy(c__ + (k - 1) * 80, ch__1, (ftnlen)8, (ftnlen)8);
	    s_copy(c__ + ((k - 1) * 80 + 3), id + (k << 3), (ftnlen)8, (
		    ftnlen)8);
	    d__1 = m[k] * k_2__;
	    mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	    s_copy(c__ + ((k - 1) * 80 + 11), ch__1, (ftnlen)8, (ftnlen)8);
	    d__1 = s[k * 3 + 1] * k_2__;
	    mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	    s_copy(c__ + ((k - 1) * 80 + 19), ch__1, (ftnlen)8, (ftnlen)8);
	    d__1 = s[k * 3 + 2] * k_2__;
	    mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	    s_copy(c__ + ((k - 1) * 80 + 27), ch__1, (ftnlen)8, (ftnlen)8);
	    d__1 = s[k * 3 + 3] * k_2__;
	    mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	    s_copy(c__ + ((k - 1) * 80 + 35), ch__1, (ftnlen)8, (ftnlen)8);
	    d__1 = rho[k] / rhocgs;
	    mio_fl2c__(ch__1, (ftnlen)8, &d__1);
	    s_copy(c__ + ((k - 1) * 80 + 43), ch__1, (ftnlen)8, (ftnlen)8);
	}

/* Write compressed output to file */
	s_wsfe(&io___928);
	do_fio(&c__1, "\f", (ftnlen)1);
	do_fio(&c__1, "6a", (ftnlen)2);
	do_fio(&c__1, (char *)&(*algor), (ftnlen)sizeof(integer));
	do_fio(&c__1, header, (ftnlen)62);
	do_fio(&c__1, (char *)&opt[4], (ftnlen)sizeof(integer));
	e_wsfe();
	i__1 = *nbod;
	for (k = 2; k <= i__1; ++k) {
	    s_wsfe(&io___929);
	    do_fio(&c__1, c__ + (k - 1) * 80, (ftnlen)51);
	    e_wsfe();
	}
    }

/* ------------------------------------------------------------------------------ */

/*  NORMAL  OUTPUT  PROCEDURE */

/* Compose a header line containing the time and number of objects */
    mio_fl2c__(ch__1, (ftnlen)8, time);
    s_copy(header, ch__1, (ftnlen)8, (ftnlen)8);
    d__1 = (doublereal) (*nbig - 1);
    mio_re2c__(ch__1, (ftnlen)8, &d__1, &c_b98, &c_b99);
    s_copy(header + 8, ch__1, (ftnlen)8, (ftnlen)8);
    d__1 = (doublereal) (*nbod - *nbig);
    mio_re2c__(ch__1, (ftnlen)8, &d__1, &c_b98, &c_b99);
    s_copy(header + 11, ch__1, (ftnlen)8, (ftnlen)8);

/* Calculate output variables for each body and convert to compressed format */
    i__1 = *nbod;
    for (k = 2; k <= i__1; ++k) {
	mco_x2ov__(rcen, rmax, &m[1], &m[k], &xh[k * 3 + 1], &xh[k * 3 + 2], &
		xh[k * 3 + 3], &vh[k * 3 + 1], &vh[k * 3 + 2], &vh[k * 3 + 3],
		 &fr, &theta, &phi, &fv, &vtheta, &vphi);

/* Object's index number and output variables */
	d__1 = (doublereal) (k - 1);
	mio_re2c__(ch__1, (ftnlen)8, &d__1, &c_b98, &c_b99);
	s_copy(c__ + (k - 1) * 80, ch__1, (ftnlen)8, (ftnlen)8);
	mio_re2c__(ch__1, (ftnlen)8, &fr, &c_b98, &rfac);
	s_copy(c__ + ((k - 1) * 80 + 3), ch__1, (ftnlen)8, (ftnlen)8);
	i__2 = nchar + 3;
	mio_re2c__(ch__1, (ftnlen)8, &theta, &c_b98, &c_b196);
	s_copy(c__ + ((k - 1) * 80 + i__2), ch__1, nchar + 11 - i__2, (ftnlen)
		8);
	i__2 = (nchar << 1) + 3;
	mio_re2c__(ch__1, (ftnlen)8, &phi, &c_b98, &c_b58);
	s_copy(c__ + ((k - 1) * 80 + i__2), ch__1, (nchar << 1) + 11 - i__2, (
		ftnlen)8);
	i__2 = nchar * 3 + 3;
	mio_re2c__(ch__1, (ftnlen)8, &fv, &c_b98, &c_b150);
	s_copy(c__ + ((k - 1) * 80 + i__2), ch__1, nchar * 3 + 11 - i__2, (
		ftnlen)8);
	i__2 = (nchar << 2) + 3;
	mio_re2c__(ch__1, (ftnlen)8, &vtheta, &c_b98, &c_b196);
	s_copy(c__ + ((k - 1) * 80 + i__2), ch__1, (nchar << 2) + 11 - i__2, (
		ftnlen)8);
	i__2 = nchar * 5 + 3;
	mio_re2c__(ch__1, (ftnlen)8, &vphi, &c_b98, &c_b58);
	s_copy(c__ + ((k - 1) * 80 + i__2), ch__1, nchar * 5 + 11 - i__2, (
		ftnlen)8);
    }

/* Write compressed output to file */
    s_wsfe(&io___936);
    do_fio(&c__1, "\f", (ftnlen)1);
    do_fio(&c__1, "6b", (ftnlen)2);
    do_fio(&c__1, header, (ftnlen)14);
    e_wsfe();
    i__1 = *nbod;
    for (k = 2; k <= i__1; ++k) {
	s_wsfe(&io___937);
	do_fio(&c__1, c__ + (k - 1) * 80, len);
	e_wsfe();
    }

    cl__1.cerr = 0;
    cl__1.cunit = 21;
    cl__1.csta = 0;
    f_clos(&cl__1);
    *opflag = 0;

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mio_out__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MIO_RE2C.FOR    (ErikSoft  27 June 1999) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Converts a REAL*8 variable X, where XMIN <= X < XMAX, into an ASCII string */
/* of 8 characters, using the new format compression: */

/* X is first converted to base 224, and then each base 224 digit is converted */
/* to an ASCII character, such that 0 -> character 32, 1 -> character 33... */
/* and 223 -> character 255. */

/* ASCII characters 0 - 31 (CTRL characters) are not used, because they */
/* cause problems when using some operating systems. */

/* ------------------------------------------------------------------------------ */

/* Character */ VOID mio_re2c__(char *ret_val, ftnlen ret_val_len, doublereal 
	*x, doublereal *xmin, doublereal *xmax)
{
    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double d_mod(doublereal *, doublereal *);

    /* Local variables */
    static integer j;
    static doublereal y, z__;



/* Input/output */

/* Local */

/* ------------------------------------------------------------------------------ */

    s_copy(ret_val, "        ", (ftnlen)8, (ftnlen)8);
    y = (*x - *xmin) / (*xmax - *xmin);

    if (y >= 1.) {
	for (j = 1; j <= 8; ++j) {
	    *(unsigned char *)&ret_val[j - 1] = 255;
	}
    } else if (y > 0.) {
	z__ = y;
	for (j = 1; j <= 8; ++j) {
	    z__ = d_mod(&z__, &c_b150) * 224.;
	    *(unsigned char *)&ret_val[j - 1] = (char) ((integer) z__ + 32);
	}
    }

/* ------------------------------------------------------------------------------ */

    return ;
} /* mio_re2c__ */


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

/*      MXX_EJEC.FOR    (ErikSoft   2 November 2000) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates the distance from the central body of each object with index */
/* I >= I0. If this distance exceeds RMAX, the object is flagged for ejection */
/* (STAT set to -3). If any object is to be ejected, EJFLAG = 1 on exit, */
/* otherwise EJFLAG = 0. */

/* Also updates the values of EN(3) and AM(3)---the change in energy and */
/* angular momentum due to collisions and ejections. */


/* N.B. All coordinates must be with respect to the central body!! */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mxx_ejec__(doublereal *time, doublereal *tstart, 
	doublereal *rmax, doublereal *en, doublereal *am, doublereal *jcen, 
	integer *i0, integer *nbod, integer *nbig, doublereal *m, doublereal *
	x, doublereal *v, doublereal *s, integer *stat, char *id, integer *
	opt, integer *ejflag, char *outfile, char *mem, integer *lmem, ftnlen 
	id_len, ftnlen outfile_len, ftnlen mem_len)
{
    /* System generated locals */
    integer i__1, i__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     f_clos(cllist *);

    /* Local variables */
    static doublereal e;
    static integer j;
    static doublereal l, r2, t1;
    static integer year;
    static doublereal rmax2;
    static integer month;
    static char flost[38];
    extern /* Subroutine */ int mxx_en__(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static char tstring[6];
    extern /* Subroutine */ int mio_jd2y__(doublereal *, integer *, integer *,
	     doublereal *);

    /* Fortran I/O blocks */
    static cilist io___953 = { 0, 23, 0, flost, 0 };
    static cilist io___955 = { 0, 23, 0, flost, 0 };




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
    --en;
    --am;
    --jcen;
    id -= 8;
    --stat;
    s -= 4;
    v -= 4;
    x -= 4;
    --m;
    --opt;
    mem -= 80;
    --lmem;

    /* Function Body */
    if (*i0 <= 0) {
	*i0 = 2;
    }
    *ejflag = 0;
    rmax2 = *rmax * *rmax;

/* Calculate initial energy and angular momentum */
    mxx_en__(&jcen[1], nbod, nbig, &m[1], &x[4], &v[4], &s[4], &e, &l);

/* Flag each object which is ejected, and set its mass to zero */
    i__1 = *nbod;
    for (j = *i0; j <= i__1; ++j) {
	r2 = x[j * 3 + 1] * x[j * 3 + 1] + x[j * 3 + 2] * x[j * 3 + 2] + x[j *
		 3 + 3] * x[j * 3 + 3];
	if (r2 > rmax2) {
	    *ejflag = 1;
	    stat[j] = -3;
	    m[j] = 0.;
	    s[j * 3 + 1] = 0.;
	    s[j * 3 + 2] = 0.;
	    s[j * 3 + 3] = 0.;

/* Write message to information file */
L20:
	    o__1.oerr = 1;
	    o__1.ounit = 23;
	    o__1.ofnmlen = 80;
	    o__1.ofnm = outfile;
	    o__1.orl = 0;
	    o__1.osta = "old";
	    o__1.oacc = "append";
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    i__2 = f_open(&o__1);
	    if (i__2 != 0) {
		goto L20;
	    }
	    if (opt[3] == 1) {
		mio_jd2y__(time, &year, &month, &t1);
		s_copy(flost, "(1x,a8,a,i10,1x,i2,1x,f8.5)", (ftnlen)38, (
			ftnlen)27);
		s_wsfe(&io___953);
		do_fio(&c__1, id + (j << 3), (ftnlen)8);
		do_fio(&c__1, mem + 5440, lmem[68]);
		do_fio(&c__1, (char *)&year, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&month, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
		e_wsfe();
	    } else {
		if (opt[3] == 3) {
		    t1 = (*time - *tstart) / 365.25;
		    s_copy(tstring, mem + 160, (ftnlen)6, (ftnlen)80);
		    s_copy(flost, "(1x,a8,a,f18.7,a)", (ftnlen)38, (ftnlen)17)
			    ;
		} else {
		    if (opt[3] == 0) {
			t1 = *time;
		    }
		    if (opt[3] == 2) {
			t1 = *time - *tstart;
		    }
		    s_copy(tstring, mem + 80, (ftnlen)6, (ftnlen)80);
		    s_copy(flost, "(1x,a8,a,f18.5,a)", (ftnlen)38, (ftnlen)17)
			    ;
		}
		s_wsfe(&io___955);
		do_fio(&c__1, id + (j << 3), (ftnlen)8);
		do_fio(&c__1, mem + 5440, lmem[68]);
		do_fio(&c__1, (char *)&t1, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, tstring, (ftnlen)6);
		e_wsfe();
	    }
	    cl__1.cerr = 0;
	    cl__1.cunit = 23;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	}
    }

/* If ejections occurred, update ELOST and LLOST */
    if (*ejflag != 0) {
	mxx_en__(&jcen[1], nbod, nbig, &m[1], &x[4], &v[4], &s[4], &en[2], &
		am[2]);
	en[3] += e - en[2];
	am[3] += l - am[2];
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mxx_ejec__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*      MXX_ELIM.FOR    (ErikSoft   13 February 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Removes any objects with STAT < 0 (i.e. those that have been flagged for */
/* removal) and reindexes all the appropriate arrays for the remaining objects. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mxx_elim__(integer *nbod, integer *nbig, doublereal *m, 
	doublereal *x, doublereal *v, doublereal *s, doublereal *rho, 
	doublereal *rceh, doublereal *rcrit, doublereal *ngf, integer *stat, 
	char *id, char *mem, integer *lmem, char *outfile, integer *nelim, 
	ftnlen id_len, ftnlen mem_len, ftnlen outfile_len)
{
    /* System generated locals */
    integer i__1, i__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *), s_wsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_wsfe(void), f_clos(cllist *);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer nbigelim, j, k, l, elim[2001];

    /* Fortran I/O blocks */
    static cilist io___961 = { 0, 23, 0, "(2a)", 0 };




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

/* Find out how many objects are to be removed */
    /* Parameter adjustments */
    id -= 8;
    --stat;
    ngf -= 5;
    --rcrit;
    --rceh;
    --rho;
    s -= 4;
    v -= 4;
    x -= 4;
    --m;
    mem -= 80;
    --lmem;

    /* Function Body */
    *nelim = 0;
    nbigelim = 0;
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	if (stat[j] < 0) {
	    ++(*nelim);
	    elim[*nelim - 1] = j;
	    if (j <= *nbig) {
		++nbigelim;
	    }
	}
    }
    elim[*nelim] = *nbod + 1;

/* Eliminate unwanted objects */
    i__1 = *nelim;
    for (k = 1; k <= i__1; ++k) {
	i__2 = elim[k] - k - 1;
	for (j = elim[k - 1] - k + 1; j <= i__2; ++j) {
	    l = j + k;
	    x[j * 3 + 1] = x[l * 3 + 1];
	    x[j * 3 + 2] = x[l * 3 + 2];
	    x[j * 3 + 3] = x[l * 3 + 3];
	    v[j * 3 + 1] = v[l * 3 + 1];
	    v[j * 3 + 2] = v[l * 3 + 2];
	    v[j * 3 + 3] = v[l * 3 + 3];
	    m[j] = m[l];
	    s[j * 3 + 1] = s[l * 3 + 1];
	    s[j * 3 + 2] = s[l * 3 + 2];
	    s[j * 3 + 3] = s[l * 3 + 3];
	    rho[j] = rho[l];
	    rceh[j] = rceh[l];
	    stat[j] = stat[l];
	    s_copy(id + (j << 3), id + (l << 3), (ftnlen)8, (ftnlen)8);
	    ngf[(j << 2) + 1] = ngf[(l << 2) + 1];
	    ngf[(j << 2) + 2] = ngf[(l << 2) + 2];
	    ngf[(j << 2) + 3] = ngf[(l << 2) + 3];
	    ngf[(j << 2) + 4] = ngf[(l << 2) + 4];
	}
    }

/* Update total number of bodies and number of Big bodies */
    *nbod -= *nelim;
    *nbig -= nbigelim;

/* If no massive bodies remain, stop the integration */
    if (*nbig < 1) {
L10:
	o__1.oerr = 1;
	o__1.ounit = 23;
	o__1.ofnmlen = 80;
	o__1.ofnm = outfile;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = "append";
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L10;
	}
	s_wsfe(&io___961);
	do_fio(&c__1, mem + 6480, lmem[81]);
	do_fio(&c__1, mem + 9920, lmem[124]);
	e_wsfe();
	cl__1.cerr = 0;
	cl__1.cunit = 23;
	cl__1.csta = 0;
	f_clos(&cl__1);
	s_stop("", (ftnlen)0);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mxx_elim__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*     MXX_EN.FOR    (ErikSoft   21 February 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates the total energy and angular-momentum for a system of objects */
/* with masses M, coordinates X, velocities V and spin angular momenta S. */

/* N.B. All coordinates and velocities must be with respect to the central */
/* ===  body. */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mxx_en__(doublereal *jcen, integer *nbod, integer *nbig, 
	doublereal *m, doublereal *xh, doublereal *vh, doublereal *s, 
	doublereal *e, doublereal *l2)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal l[3], v[6000]	/* was [3][2000] */, x[6000]	/* 
	    was [3][2000] */, r2, u2, u4, u6, ke, pe, dx, dy, dz, r_1__, 
	    r_2__, r_4__, r_6__, tmp, tmp2[8000]	/* was [4][2000] */, 
	    temp;
    static integer itmp[8], iflag;
    extern /* Subroutine */ int mco_h2b__(doublereal *, doublereal *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *);



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
    --jcen;
    s -= 4;
    vh -= 4;
    xh -= 4;
    --m;

    /* Function Body */
    ke = 0.;
    pe = 0.;
    l[0] = 0.;
    l[1] = 0.;
    l[2] = 0.;

/* Convert to barycentric coordinates and velocities */
    mco_h2b__(&temp, &jcen[1], nbod, nbig, &temp, &m[1], &xh[4], &vh[4], x, v,
	     tmp2, &iflag, itmp);

/* Do the spin angular momenta first (probably the smallest terms) */
    i__1 = *nbod;
    for (j = 1; j <= i__1; ++j) {
	l[0] += s[j * 3 + 1];
	l[1] += s[j * 3 + 2];
	l[2] += s[j * 3 + 3];
    }

/* Orbital angular momentum and kinetic energy terms */
    i__1 = *nbod;
    for (j = 1; j <= i__1; ++j) {
	l[0] += m[j] * (x[j * 3 - 2] * v[j * 3 - 1] - x[j * 3 - 1] * v[j * 3 
		- 2]);
	l[1] += m[j] * (x[j * 3 - 1] * v[j * 3 - 3] - x[j * 3 - 3] * v[j * 3 
		- 1]);
	l[2] += m[j] * (x[j * 3 - 3] * v[j * 3 - 2] - x[j * 3 - 2] * v[j * 3 
		- 3]);
	ke += m[j] * (v[j * 3 - 3] * v[j * 3 - 3] + v[j * 3 - 2] * v[j * 3 - 
		2] + v[j * 3 - 1] * v[j * 3 - 1]);
    }

/* Potential energy terms due to pairs of bodies */
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	tmp = 0.;
	i__2 = *nbod;
	for (k = j + 1; k <= i__2; ++k) {
	    dx = x[k * 3 - 3] - x[j * 3 - 3];
	    dy = x[k * 3 - 2] - x[j * 3 - 2];
	    dz = x[k * 3 - 1] - x[j * 3 - 1];
	    r2 = dx * dx + dy * dy + dz * dz;
	    if (r2 != 0.) {
		tmp += m[k] / sqrt(r2);
	    }
	}
	pe -= tmp * m[j];
    }

/* Potential energy terms involving the central body */
    i__1 = *nbod;
    for (j = 2; j <= i__1; ++j) {
	dx = x[j * 3 - 3] - x[0];
	dy = x[j * 3 - 2] - x[1];
	dz = x[j * 3 - 1] - x[2];
	r2 = dx * dx + dy * dy + dz * dz;
	if (r2 != 0.) {
	    pe -= m[1] * m[j] / sqrt(r2);
	}
    }

/* Corrections for oblateness */
    if (jcen[1] != 0. || jcen[2] != 0. || jcen[3] != 0.) {
	i__1 = *nbod;
	for (j = 2; j <= i__1; ++j) {
	    r2 = xh[j * 3 + 1] * xh[j * 3 + 1] + xh[j * 3 + 2] * xh[j * 3 + 2]
		     + xh[j * 3 + 3] * xh[j * 3 + 3];
	    r_1__ = 1. / sqrt(r2);
	    r_2__ = r_1__ * r_1__;
	    r_4__ = r_2__ * r_2__;
	    r_6__ = r_4__ * r_2__;
	    u2 = xh[j * 3 + 3] * xh[j * 3 + 3] * r_2__;
	    u4 = u2 * u2;
	    u6 = u4 * u2;
	    pe += m[1] * m[j] * r_1__ * (jcen[1] * r_2__ * (u2 * 1.5 - .5) + 
		    jcen[2] * r_4__ * (u4 * 4.375 - u2 * 3.75 + .375) + jcen[
		    3] * r_6__ * (u6 * 14.4375 - u4 * 19.6875 + u2 * 6.5625 - 
		    .3125));
	}
    }

    *e = ke * .5 + pe;
    *l2 = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mxx_en__ */


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/*     MXX_JAC.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Calculates the Jacobi constant for massless particles. This assumes that */
/* there are only 2 massive bodies (including the central body) moving on */
/* circular orbits. */

/* N.B. All coordinates and velocities must be heliocentric!! */
/* === */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mxx_jac__(doublereal *jcen, integer *nbod, integer *nbig,
	 doublereal *m, doublereal *xh, doublereal *vh, doublereal *jac)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal d__;
    static integer j;
    static doublereal n, r__, v[6000]	/* was [3][2000] */, x[6000]	/* 
	    was [3][2000] */, a2, dx, dy, dz, tmp2[8000]	/* was [4][
	    2000] */, temp;
    static integer itmp[8], iflag;
    extern /* Subroutine */ int mco_h2b__(doublereal *, doublereal *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *);



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
    --jcen;
    vh -= 4;
    xh -= 4;
    --m;
    --jac;

    /* Function Body */
    mco_h2b__(&temp, &jcen[1], nbod, nbig, &temp, &m[1], &xh[4], &vh[4], x, v,
	     tmp2, &iflag, itmp);
    dx = x[3] - x[0];
    dy = x[4] - x[1];
    dz = x[5] - x[2];
    a2 = dx * dx + dy * dy + dz * dz;
    n = sqrt((m[1] + m[2]) / (a2 * sqrt(a2)));

    i__1 = *nbod;
    for (j = *nbig + 1; j <= i__1; ++j) {
	dx = x[j * 3 - 3] - x[0];
	dy = x[j * 3 - 2] - x[1];
	dz = x[j * 3 - 1] - x[2];
	r__ = sqrt(dx * dx + dy * dy + dz * dz);
	dx = x[j * 3 - 3] - x[3];
	dy = x[j * 3 - 2] - x[4];
	dz = x[j * 3 - 1] - x[5];
	d__ = sqrt(dx * dx + dy * dy + dz * dz);

	jac[j] = (v[j * 3 - 3] * v[j * 3 - 3] + v[j * 3 - 2] * v[j * 3 - 2] + 
		v[j * 3 - 1] * v[j * 3 - 1]) * .5 - m[1] / r__ - m[2] / d__ - 
		n * (x[j * 3 - 3] * v[j * 3 - 2] - x[j * 3 - 2] * v[j * 3 - 3]
		);
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mxx_jac__ */


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

/*      MXX_SYNC.FOR    (ErikSoft   2 March 2001) */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* Author: John E. Chambers */

/* Synchronizes the epochs of NBIG Big bodies (having a common epoch) and */
/* NBOD-NBIG Small bodies (possibly having differing epochs), for an */
/* integration using MERCURY. */
/* The Small bodies are picked up in order starting with the one with epoch */
/* furthest from the time, TSTART, at which the main integration will begin */
/* producing output. */

/* N.B. The synchronization integrations use Everhart's RA15 routine. */
/* --- */

/* ------------------------------------------------------------------------------ */

/* Subroutine */ int mxx_sync__(doublereal *time, doublereal *tstart, 
	doublereal *h0, doublereal *tol, doublereal *jcen, integer *nbod, 
	integer *nbig, doublereal *m, doublereal *x, doublereal *v, 
	doublereal *s, doublereal *rho, doublereal *rceh, integer *stat, char 
	*id, doublereal *epoch, doublereal *ngf, integer *opt, integer *
	ngflag, ftnlen id_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    extern /* Subroutine */ int mxx_sort__(integer *, doublereal *, integer *)
	    ;
    static doublereal h__;
    static integer j, k, l, ice[2000], jce[2000], nce;
    static doublereal hdid;
    static integer indx[2000];
    static doublereal temp;
    static integer nsml;
    static char ctemp[8*2000];
    static integer itemp, jtemp[2000];
    static doublereal epsml[2000], rcrit[2000], rtemp[2000], rphys[2000];
    static integer raflag, nsofar;
    static doublereal tsmall;
    extern /* Subroutine */ int mfo_all__(doublereal *, doublereal *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, integer *, integer *), 
	    mdt_ra15__(doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, S_fp);



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

/* Reorder Small bodies by epoch so that ep(1) is furthest from TSTART */
    /* Parameter adjustments */
    --jcen;
    ngf -= 5;
    --epoch;
    id -= 8;
    --stat;
    --rceh;
    --rho;
    s -= 4;
    v -= 4;
    x -= 4;
    --m;
    --opt;

    /* Function Body */
    nsml = *nbod - *nbig;
    i__1 = *nbod;
    for (j = *nbig + 1; j <= i__1; ++j) {
	epsml[j - *nbig - 1] = epoch[j];
    }
    mxx_sort__(&nsml, epsml, indx);

    if ((d__1 = epsml[0] - *tstart, abs(d__1)) < (d__2 = epsml[nsml - 1] - *
	    tstart, abs(d__2))) {
	k = nsml + 1;
	i__1 = nsml / 2;
	for (j = 1; j <= i__1; ++j) {
	    l = k - j;
	    temp = epsml[j - 1];
	    epsml[j - 1] = epsml[l - 1];
	    epsml[l - 1] = temp;
	    itemp = indx[j - 1];
	    indx[j - 1] = indx[l - 1];
	    indx[l - 1] = itemp;
	}
    }

    i__1 = *nbod;
    for (j = *nbig + 1; j <= i__1; ++j) {
	epoch[j] = epsml[j - *nbig - 1];
    }

/* Reorder the other arrays associated with each Small body */
    for (k = 1; k <= 3; ++k) {
	i__1 = nsml;
	for (j = 1; j <= i__1; ++j) {
	    rtemp[j - 1] = x[k + (j + *nbig) * 3];
	}
	i__1 = nsml;
	for (j = 1; j <= i__1; ++j) {
	    x[k + (j + *nbig) * 3] = rtemp[indx[j - 1] - 1];
	}
	i__1 = nsml;
	for (j = 1; j <= i__1; ++j) {
	    rtemp[j - 1] = v[k + (j + *nbig) * 3];
	}
	i__1 = nsml;
	for (j = 1; j <= i__1; ++j) {
	    v[k + (j + *nbig) * 3] = rtemp[indx[j - 1] - 1];
	}
	i__1 = nsml;
	for (j = 1; j <= i__1; ++j) {
	    rtemp[j - 1] = s[k + (j + *nbig) * 3];
	}
	i__1 = nsml;
	for (j = 1; j <= i__1; ++j) {
	    s[k + (j + *nbig) * 3] = rtemp[indx[j - 1] - 1];
	}
    }

    i__1 = nsml;
    for (j = 1; j <= i__1; ++j) {
	rtemp[j - 1] = m[j + *nbig];
    }
    i__1 = nsml;
    for (j = 1; j <= i__1; ++j) {
	m[j + *nbig] = rtemp[indx[j - 1] - 1];
    }
    i__1 = nsml;
    for (j = 1; j <= i__1; ++j) {
	rtemp[j - 1] = rceh[j + *nbig];
    }
    i__1 = nsml;
    for (j = 1; j <= i__1; ++j) {
	rceh[j + *nbig] = rtemp[indx[j - 1] - 1];
    }
    i__1 = nsml;
    for (j = 1; j <= i__1; ++j) {
	rtemp[j - 1] = rho[j + *nbig];
    }
    i__1 = nsml;
    for (j = 1; j <= i__1; ++j) {
	rho[j + *nbig] = rtemp[indx[j - 1] - 1];
    }

    i__1 = nsml;
    for (j = 1; j <= i__1; ++j) {
	s_copy(ctemp + (j - 1 << 3), id + (j + *nbig << 3), (ftnlen)8, (
		ftnlen)8);
	jtemp[j - 1] = stat[j + *nbig];
    }
    i__1 = nsml;
    for (j = 1; j <= i__1; ++j) {
	s_copy(id + (j + *nbig << 3), ctemp + (indx[j - 1] - 1 << 3), (ftnlen)
		8, (ftnlen)8);
	stat[j + *nbig] = jtemp[indx[j - 1] - 1];
    }

/* Integrate Small bodies up to the same epoch */
    h__ = *h0;
    tsmall = *h0 * 1e-12;
    raflag = 0;

    i__1 = *nbod;
    for (j = *nbig + 1; j <= i__1; ++j) {
	nsofar = j - 1;
	while((d__1 = *time - epoch[j], abs(d__1)) > tsmall) {
	    temp = epoch[j] - *time;
/* Computing MAX */
/* Computing MIN */
	    d__3 = abs(temp), d__4 = abs(h__);
	    d__2 = min(d__3,d__4);
	    d__1 = max(d__2,tsmall);
	    h__ = d_sign(&d__1, &temp);
	    mdt_ra15__(time, &h__, &hdid, tol, &jcen[1], &nsofar, nbig, &m[1],
		     &x[4], &v[4], &s[4], rphys, rcrit, &ngf[5], &stat[1], &
		    raflag, ngflag, &opt[1], &nce, ice, jce, (S_fp)mfo_all__);
	    *time += hdid;
	}
	raflag = 1;
    }

/* ------------------------------------------------------------------------------ */

    return 0;
} /* mxx_sync__ */


/* ************************************************************************* */
/*                        DRIFT_DAN.F */
/* ************************************************************************* */
/* This subroutine does the Danby and decides which vbles to use */

/*             Input: */
/*                 nbod          ==>  number of massive bodies (int scalar) */
/*                 mass          ==>  mass of bodies (real array) */
/*                 x0,y0,z0         ==>  initial position in jacobi coord */
/*                                    (real scalar) */
/*                 vx0,vy0,vz0      ==>  initial position in jacobi coord */
/*                                    (real scalar) */
/*                 dt0            ==>  time step */
/*             Output: */
/*                 x0,y0,z0         ==>  final position in jacobi coord */
/*                                       (real scalars) */
/*                 vx0,vy0,vz0      ==>  final position in jacobi coord */
/*                                       (real scalars) */
/*                 iflg             ==>  integer flag (zero if satisfactory) */
/* 					      (non-zero if nonconvergence) */

/* Authors:  Hal Levison & Martin Duncan */
/* Date:    2/10/93 */
/* Last revision: April 6/93 - MD adds dt and keeps dt0 unchanged */
/* Subroutine */ int drift_dan__(doublereal *mu, doublereal *x0, doublereal *
	y0, doublereal *z0, doublereal *vx0, doublereal *vy0, doublereal *vz0,
	 doublereal *dt0, integer *iflg)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal a, c__, f, g, s, u, x, y, z__, c1, c2, c3, r0, ec, dm, 
	    en, fp, dt, es, vx, vy, vz;
    extern /* Subroutine */ int drift_kepu__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static doublereal v0s, asq, esq;
    extern /* Subroutine */ int drift_kepmd__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal fchk, fdot, gdot, xkep, alpha;

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
/* ...  Inputs and Outputs: */
/* ...  Output */
/* ...  Internals: */
/* ---- */
/* ...  Executable code */
/* ...  Set dt = dt0 to be sure timestep is not altered while solving */
/* ...  for new coords. */
    dt = *dt0;
    *iflg = 0;
    r0 = sqrt(*x0 * *x0 + *y0 * *y0 + *z0 * *z0);
    v0s = *vx0 * *vx0 + *vy0 * *vy0 + *vz0 * *vz0;
    u = *x0 * *vx0 + *y0 * *vy0 + *z0 * *vz0;
    alpha = *mu * 2.f / r0 - v0s;
    if (alpha > 0.) {
	a = *mu / alpha;
	asq = a * a;
	en = sqrt(*mu / (a * asq));
	ec = 1. - r0 / a;
	es = u / (en * asq);
	esq = ec * ec + es * es;
	dm = dt * en - (integer) (dt * en / 6.2831853071795862) * 
		6.2831853071795862;
	dt = dm / en;
	if (dm * dm > .16 || esq > .36) {
	    goto L100;
	}
	if (esq * dm * dm < .0016f) {
	    drift_kepmd__(&dm, &es, &ec, &xkep, &s, &c__);
	    fchk = xkep - ec * s + es * (1.f - c__) - dm;
	    if (fchk * fchk > 1e-13) {
		*iflg = 1;
		return 0;
	    }
	    fp = 1.f - ec * c__ + es * s;
	    f = a / r0 * (c__ - 1.f) + 1.f;
	    g = dt + (s - xkep) / en;
	    fdot = -(a / (r0 * fp)) * en * s;
	    gdot = (c__ - 1.f) / fp + 1.f;
	    x = *x0 * f + *vx0 * g;
	    y = *y0 * f + *vy0 * g;
	    z__ = *z0 * f + *vz0 * g;
	    vx = *x0 * fdot + *vx0 * gdot;
	    vy = *y0 * fdot + *vy0 * gdot;
	    vz = *z0 * fdot + *vz0 * gdot;
	    *x0 = x;
	    *y0 = y;
	    *z0 = z__;
	    *vx0 = vx;
	    *vy0 = vy;
	    *vz0 = vz;
	    *iflg = 0;
	    return 0;
	}
    }
L100:
    drift_kepu__(&dt, &r0, mu, &alpha, &u, &fp, &c1, &c2, &c3, iflg);
    if (*iflg == 0) {
	f = 1.f - *mu / r0 * c2;
	g = dt - *mu * c3;
	fdot = -(*mu / (fp * r0)) * c1;
	gdot = 1.f - *mu / fp * c2;
	x = *x0 * f + *vx0 * g;
	y = *y0 * f + *vy0 * g;
	z__ = *z0 * f + *vz0 * g;
	vx = *x0 * fdot + *vx0 * gdot;
	vy = *y0 * fdot + *vy0 * gdot;
	vz = *z0 * fdot + *vz0 * gdot;
	*x0 = x;
	*y0 = y;
	*z0 = z__;
	*vx0 = vx;
	*vy0 = vy;
	*vz0 = vz;
    }
    return 0;
} /* drift_dan__ */


/* ********************************************************************# */
/*                  DRIFT_KEPMD */
/* ********************************************************************# */
/*  Subroutine for solving kepler's equation in difference form for an */
/*  ellipse, given SMALL dm and SMALL eccentricity.  See DRIFT_DAN.F */
/*  for the criteria. */
/*  WARNING - BUILT FOR SPEED : DOES NOT CHECK HOW WELL THE ORIGINAL */
/*  EQUATION IS SOLVED! (CAN DO THAT IN THE CALLING ROUTINE BY */
/*  CHECKING HOW CLOSE (x - ec*s +es*(1.-c) - dm) IS TO ZERO. */

/* 	Input: */
/* 	    dm		==> increment in mean anomaly M (real*8 scalar) */
/* 	    es,ec       ==> ecc. times sin and cos of E_0 (real*8 scalars) */

/*       Output: */
/*            x          ==> solution to Kepler's difference eqn (real*8 scalar) */
/*            s,c        ==> sin and cosine of x (real*8 scalars) */

/* drift_dan */
/* Subroutine */ int drift_kepmd__(doublereal *dm, doublereal *es, doublereal 
	*ec, doublereal *x, doublereal *s, doublereal *c__)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal f, q, y, fp, dx, fpp, fac1, fac2, fppp;

/* ...    Inputs */
/* ...	Outputs */
/* ...    Internals */
/* ...    calc initial guess for root */
    fac1 = 1. / (1. - *ec);
    q = fac1 * *dm;
    fac2 = *es * *es * fac1 - *ec / 3.;
    *x = q * (1. - fac1 * .5 * q * (*es - q * fac2));
/* ...  excellent approx. to sin and cos of x for small x. */
    y = *x * *x;
    *s = *x * (39916800. - y * (6652800. - y * (332640. - y * (7920. - y * (
	    110. - y))))) / 39916800.;
    *c__ = sqrt(1. - *s * *s);
/* ...    Compute better value for the root using quartic Newton method */
    f = *x - *ec * *s + *es * (1.f - *c__) - *dm;
    fp = 1.f - *ec * *c__ + *es * *s;
    fpp = *ec * *s + *es * *c__;
    fppp = *ec * *c__ - *es * *s;
    dx = -f / fp;
    dx = -f / (fp + dx * .5f * fpp);
    dx = -f / (fp + dx * .5f * fpp + dx * .16666666666666666f * dx * fppp);
    *x += dx;
/* ...  excellent approx. to sin and cos of x for small x. */
    y = *x * *x;
    *s = *x * (39916800. - y * (6652800. - y * (332640. - y * (7920. - y * (
	    110. - y))))) / 39916800.;
    *c__ = sqrt(1. - *s * *s);
    return 0;
} /* drift_kepmd__ */

/* ----------------------------------------------------------------------------- */

/* ************************************************************************* */
/*                        DRIFT_KEPU.F */
/* ************************************************************************* */
/* subroutine for solving kepler's equation using universal variables. */

/*             Input: */
/*                 dt            ==>  time step (real scalor) */
/*                 r0            ==>  Distance between `Sun' and paritcle */
/*                                     (real scalor) */
/*                 mu            ==>  Reduced mass of system (real scalor) */
/*                 alpha         ==>  energy (real scalor) */
/*                 u             ==>  angular momentun  (real scalor) */
/*             Output: */
/*                 fp            ==>  f' from p170 */
/*                                       (real scalors) */
/*                 c1,c2,c3      ==>  c's from p171-172 */
/*                                       (real scalors) */
/*                 iflg          ==>  =0 if converged; !=0 if not */

/* Author:  Hal Levison */
/* Date:    2/3/93 */
/* Last revision: 2/3/93 */
/* Subroutine */ int drift_kepu__(doublereal *dt, doublereal *r0, doublereal *
	mu, doublereal *alpha, doublereal *u, doublereal *fp, doublereal *c1, 
	doublereal *c2, doublereal *c3, integer *iflg)
{
    extern /* Subroutine */ int drift_kepu_guess__(doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal s, fn, fo, st;
    extern /* Subroutine */ int drift_kepu_lag__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , drift_kepu_new__(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *), 
	    drift_kepu_fchk__(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);

/* ...  Inputs: */
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
/* ...  Outputs: */
/* ...  Internals: */
/* ---- */
/* ...  Executable code */
    drift_kepu_guess__(dt, r0, mu, alpha, u, &s);
    st = s;
/* ..     store initial guess for possible use later in */
/* ..     laguerre's method, in case newton's method fails. */
    drift_kepu_new__(&s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg);
    if (*iflg != 0) {
	drift_kepu_fchk__(dt, r0, mu, alpha, u, &st, &fo);
	drift_kepu_fchk__(dt, r0, mu, alpha, u, &s, &fn);
	if (abs(fo) < abs(fn)) {
	    s = st;
	}
	drift_kepu_lag__(&s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg);
    }
    return 0;
} /* drift_kepu__ */

/* ---------------------------------------------------------------------- */

/* ************************************************************************* */
/*                        DRIFT_KEPU_FCHK.F */
/* ************************************************************************* */
/* Returns the value of the function f of which we are trying to find the root */
/* in universal variables. */

/*             Input: */
/*                 dt            ==>  time step (real scalar) */
/*                 r0            ==>  Distance between `Sun' and particle */
/*                                     (real scalar) */
/*                 mu            ==>  Reduced mass of system (real scalar) */
/*                 alpha         ==>  Twice the binding energy (real scalar) */
/*                 u             ==>  Vel. dot radial vector (real scalar) */
/*                 s             ==>  Approx. root of f */
/*             Output: */
/*                 f             ==>  function value ( = 0 if O.K.) (integer) */

/* Author:  Martin Duncan */
/* Date:    March 12/93 */
/* Last revision: March 12/93 */
/* drift_kepu */
/* Subroutine */ int drift_kepu_fchk__(doublereal *dt, doublereal *r0, 
	doublereal *mu, doublereal *alpha, doublereal *u, doublereal *s, 
	doublereal *f)
{
    static doublereal x, c0, c1, c2, c3;
    extern /* Subroutine */ int drift_kepu_stumpff__(doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *);

/* ...  Inputs: */
/* ...  Outputs: */
/* ...  Internals: */
/* ---- */
/* ...  Executable code */
    x = *s * *s * *alpha;
    drift_kepu_stumpff__(&x, &c0, &c1, &c2, &c3);
    c1 *= *s;
    c2 = c2 * *s * *s;
    c3 = c3 * *s * *s * *s;
    *f = *r0 * c1 + *u * c2 + *mu * c3 - *dt;
    return 0;
} /* drift_kepu_fchk__ */

/* ------------------------------------------------------------------- */

/* ************************************************************************* */
/*                        DRIFT_KEPU_GUESS.F */
/* ************************************************************************* */
/* Initial guess for solving kepler's equation using universal variables. */

/*             Input: */
/*                 dt            ==>  time step (real scalor) */
/*                 r0            ==>  Distance between `Sun' and paritcle */
/*                                     (real scalor) */
/*                 mu            ==>  Reduced mass of system (real scalor) */
/*                 alpha         ==>  energy (real scalor) */
/*                 u             ==>  angular momentun  (real scalor) */
/*             Output: */
/*                 s             ==>  initial guess for the value of */
/*                                    universal variable */

/* Author:  Hal Levison & Martin Duncan */
/* Date:    3/12/93 */
/* Last revision: April 6/93 */
/* Modified by JEC: 8/6/98 */
/*   drift_kepu_fchk */
/* Subroutine */ int drift_kepu_guess__(doublereal *dt, doublereal *r0, 
	doublereal *mu, doublereal *alpha, doublereal *u, doublereal *s)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    extern /* Subroutine */ int mco_sine__(doublereal *, doublereal *, 
	    doublereal *);
    static doublereal a, e, x, y, ec, en, es, cy, sy;
    extern /* Subroutine */ int drift_kepu_p3solve__(doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *);
    static integer iflg;
    static doublereal sigma;

/* ...  Inputs: */
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
/* ...  Inputs and Outputs: */
/* ...  Internals: */
/* ---- */
/* ...  Executable code */
    if (*alpha > 0.f) {
/* ...       find initial guess for elliptic motion */
	if (*dt / *r0 <= .4f) {
	    *s = *dt / *r0 - *dt * *dt * *u / (*r0 * 2.f * *r0 * *r0);
	    return 0;
	} else {
	    a = *mu / *alpha;
	    en = sqrt(*mu / (a * a * a));
	    ec = 1.f - *r0 / a;
	    es = *u / (en * a * a);
	    e = sqrt(ec * ec + es * es);
	    y = en * *dt - es;

	    mco_sine__(&y, &sy, &cy);

	    d__1 = es * cy + ec * sy;
	    sigma = d_sign(&c_b150, &d__1);
	    x = y + sigma * .85f * e;
	    *s = x / sqrt(*alpha);
	}
    } else {
/* ...       find initial guess for hyperbolic motion. */
	drift_kepu_p3solve__(dt, r0, mu, alpha, u, s, &iflg);
	if (iflg != 0) {
	    *s = *dt / *r0;
	}
    }
    return 0;
} /* drift_kepu_guess__ */

/* ------------------------------------------------------------------- */

/* ************************************************************************* */
/*                        DRIFT_KEPU_LAG.F */
/* ************************************************************************* */
/* subroutine for solving kepler's equation in universal variables. */
/* using LAGUERRE'S METHOD */

/*             Input: */
/*                 s             ==>  inital value of universal variable */
/*                 dt            ==>  time step (real scalor) */
/*                 r0            ==>  Distance between `Sun' and paritcle */
/*                                     (real scalor) */
/*                 mu            ==>  Reduced mass of system (real scalor) */
/*                 alpha         ==>  energy (real scalor) */
/*                 u             ==>  angular momentun  (real scalor) */
/*             Output: */
/*                 s             ==>  final value of universal variable */
/*                 fp            ==>  f' from p170 */
/*                                       (real scalors) */
/*                 c1,c2,c3      ==>  c's from p171-172 */
/*                                       (real scalors) */
/*                 iflgn          ==>  =0 if converged; !=0 if not */

/* Author:  Hal Levison */
/* Date:    2/3/93 */
/* Last revision: 4/21/93 */
/*   drift_kepu_guess */
/* Subroutine */ int drift_kepu_lag__(doublereal *s, doublereal *dt, 
	doublereal *r0, doublereal *mu, doublereal *alpha, doublereal *u, 
	doublereal *fp, doublereal *c1, doublereal *c2, doublereal *c3, 
	integer *iflg)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal f, x, c0;
    static integer nc;
    static doublereal ds, ln, fdt, fpp;
    extern /* Subroutine */ int drift_kepu_stumpff__(doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *);
    static integer ncmax;

/* ...  Inputs: */
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
/* ...  Outputs: */
/* ...  Internals: */
/* ---- */
/* ...  Executable code */
/* ...    To get close approch needed to take lots of iterations if alpha<0 */
    if (*alpha < 0.f) {
	ncmax = 400;
    } else {
	ncmax = 400;
    }
    ln = 5.f;
/* ...    start laguere's method */
    i__1 = ncmax;
    for (nc = 0; nc <= i__1; ++nc) {
	x = *s * *s * *alpha;
	drift_kepu_stumpff__(&x, &c0, c1, c2, c3);
	*c1 *= *s;
	*c2 = *c2 * *s * *s;
	*c3 = *c3 * *s * *s * *s;
	f = *r0 * *c1 + *u * *c2 + *mu * *c3 - *dt;
	*fp = *r0 * c0 + *u * *c1 + *mu * *c2;
	fpp = (*alpha * -40.f + *mu) * *c1 + *u * c0;
	ds = -ln * f / (*fp + d_sign(&c_b150, fp) * sqrt((d__1 = (ln - 1.f) * 
		(ln - 1.f) * *fp * *fp - (ln - 1.f) * ln * f * fpp, abs(d__1))
		));
	*s += ds;
	fdt = f / *dt;
/* ..        quartic convergence */
	if (fdt * fdt < 1e-26) {
	    *iflg = 0;
	    return 0;
	}
/* ...      Laguerre's method succeeded */
    }
    *iflg = 2;
    return 0;
} /* drift_kepu_lag__ */

/* ----------------------------------------------------------------------- */

/* ************************************************************************* */
/*                        DRIFT_KEPU_NEW.F */
/* ************************************************************************* */
/* subroutine for solving kepler's equation in universal variables. */
/* using NEWTON'S METHOD */

/*             Input: */
/*                 s             ==>  inital value of universal variable */
/*                 dt            ==>  time step (real scalor) */
/*                 r0            ==>  Distance between `Sun' and paritcle */
/*                                     (real scalor) */
/*                 mu            ==>  Reduced mass of system (real scalor) */
/*                 alpha         ==>  energy (real scalor) */
/*                 u             ==>  angular momentun  (real scalor) */
/*             Output: */
/*                 s             ==>  final value of universal variable */
/*                 fp            ==>  f' from p170 */
/*                                       (real scalors) */
/*                 c1,c2,c3      ==>  c's from p171-172 */
/*                                       (real scalors) */
/*                 iflgn          ==>  =0 if converged; !=0 if not */

/* Author:  Hal Levison */
/* Date:    2/3/93 */
/* Last revision: 4/21/93 */
/* Modified by JEC: 31/3/98 */
/*    drift_kepu_leg */
/* Subroutine */ int drift_kepu_new__(doublereal *s, doublereal *dt, 
	doublereal *r0, doublereal *mu, doublereal *alpha, doublereal *u, 
	doublereal *fp, doublereal *c1, doublereal *c2, doublereal *c3, 
	integer *iflgn)
{
    static doublereal f, x, c0, s2;
    static integer nc;
    static doublereal ds, fdt, fpp;
    extern /* Subroutine */ int drift_kepu_stumpff__(doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *);
    static doublereal fppp;

/* ...  Inputs: */
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
/* ...  Outputs: */
/* ...  Internals: */
/* ---- */
/* ...  Executable code */
    for (nc = 0; nc <= 6; ++nc) {
	s2 = *s * *s;
	x = s2 * *alpha;
	drift_kepu_stumpff__(&x, &c0, c1, c2, c3);
	*c1 *= *s;
	*c2 *= s2;
	*c3 = *c3 * *s * s2;
	f = *r0 * *c1 + *u * *c2 + *mu * *c3 - *dt;
	*fp = *r0 * c0 + *u * *c1 + *mu * *c2;
	fpp = (*mu - *r0 * *alpha) * *c1 + *u * c0;
	fppp = (*mu - *r0 * *alpha) * c0 - *u * *alpha * *c1;
	ds = -f / *fp;
	ds = -f / (*fp + ds * .5 * fpp);
	ds = -f / (*fp + ds * .5 * fpp + ds * ds * fppp * .1666666666666667);
	*s += ds;
	fdt = f / *dt;
/* ..      quartic convergence */
	if (fdt * fdt < 1e-26) {
	    *iflgn = 0;
	    return 0;
	}
/* ...     newton's method succeeded */
    }
/* ..     newton's method failed */
    *iflgn = 1;
    return 0;
} /* drift_kepu_new__ */

/* ---------------------------------------------------------------------- */

/* ************************************************************************* */
/*                        DRIFT_KEPU_P3SOLVE.F */
/* ************************************************************************* */
/* Returns the real root of cubic often found in solving kepler */
/* problem in universal variables. */

/*             Input: */
/*                 dt            ==>  time step (real scalar) */
/*                 r0            ==>  Distance between `Sun' and paritcle */
/*                                     (real scalar) */
/*                 mu            ==>  Reduced mass of system (real scalar) */
/*                 alpha         ==>  Twice the binding energy (real scalar) */
/*                 u             ==>  Vel. dot radial vector (real scalar) */
/*             Output: */
/*                 s             ==>  solution of cubic eqn for the */
/*                                    universal variable */
/*                 iflg          ==>  success flag ( = 0 if O.K.) (integer) */

/* Author:  Martin Duncan */
/* Date:    March 12/93 */
/* Last revision: March 12/93 */
/* drift_kepu_new */
/* Subroutine */ int drift_kepu_p3solve__(doublereal *dt, doublereal *r0, 
	doublereal *mu, doublereal *alpha, doublereal *u, doublereal *s, 
	integer *iflg)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal q, r__, a0, a1, a2, p1, p2, sq, sq2, denom;

/* ...  Inputs: */
/* ...  Outputs: */
/* ...  Internals: */
/* ---- */
/* ...  Executable code */
    denom = (*mu - *alpha * *r0) / 6.;
    a2 = *u * .5f / denom;
    a1 = *r0 / denom;
    a0 = -(*dt) / denom;
    q = (a1 - a2 * a2 / 3.) / 3.;
/* Computing 3rd power */
    d__1 = a2;
    r__ = (a1 * a2 - a0 * 3.) / 6. - d__1 * (d__1 * d__1) / 27.;
/* Computing 3rd power */
    d__1 = q;
/* Computing 2nd power */
    d__2 = r__;
    sq2 = d__1 * (d__1 * d__1) + d__2 * d__2;
    if (sq2 >= 0.) {
	sq = sqrt(sq2);
	if (r__ + sq <= 0.) {
	    d__1 = -(r__ + sq);
	    p1 = -pow_dd(&d__1, &c_b1082);
	} else {
	    d__1 = r__ + sq;
	    p1 = pow_dd(&d__1, &c_b1082);
	}
	if (r__ - sq <= 0.) {
	    d__1 = -(r__ - sq);
	    p2 = -pow_dd(&d__1, &c_b1082);
	} else {
	    d__1 = r__ - sq;
	    p2 = pow_dd(&d__1, &c_b1082);
	}
	*iflg = 0;
	*s = p1 + p2 - a2 / 3.;
    } else {
	*iflg = 1;
	*s = 0.;
    }
    return 0;
} /* drift_kepu_p3solve__ */

/* ------------------------------------------------------------------- */

/* ************************************************************************* */
/*                        DRIFT_KEPU_STUMPFF.F */
/* ************************************************************************* */
/* subroutine for the calculation of stumpff functions */
/* see Danby p.172  equations 6.9.15 */

/*             Input: */
/*                 x             ==>  argument */
/*             Output: */
/*                 c0,c1,c2,c3   ==>  c's from p171-172 */
/*                                       (real scalors) */
/* Author:  Hal Levison */
/* Date:    2/3/93 */
/* Last revision: 2/3/93 */
/* Modified by JEC: 31/3/98 */

/*   drift_kepu_p3solve */
/* Subroutine */ int drift_kepu_stumpff__(doublereal *x, doublereal *c0, 
	doublereal *c1, doublereal *c2, doublereal *c3)
{
    static integer i__, n;
    static doublereal x2, x3, x4, x5, x6, xm;

/* ...  Inputs: */
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
/* ...  Outputs: */
/* ...  Internals: */
/* ---- */
/* ...  Executable code */
    n = 0;
    xm = .1f;
    while(abs(*x) >= xm) {
	++n;
	*x *= .25;
    }

    x2 = *x * *x;
    x3 = *x * x2;
    x4 = x2 * x2;
    x5 = x2 * x3;
    x6 = x3 * x3;

    *c2 = x6 * 1.147074559772972e-11 - x5 * 2.08767569878681e-9 + x4 * 
	    2.755731922398589e-7 - x3 * 2.48015873015873e-5 + x2 * 
	    .001388888888888889 - *x * .04166666666666667 + .5;

    *c3 = x6 * 7.647163731819816e-13 - x5 * 1.605904383682161e-10 + x4 * 
	    2.505210838544172e-8 - x3 * 2.755731922398589e-6 + x2 * 
	    1.984126984126984e-4 - *x * .008333333333333333 + 
	    .1666666666666667;

    *c1 = 1.f - *x * *c3;
    *c0 = 1.f - *x * *c2;

    if (n != 0) {
	for (i__ = n; i__ >= 1; --i__) {
	    *c3 = (*c2 + *c0 * *c3) * .25;
	    *c2 = *c1 * *c1 * .5;
	    *c1 = *c0 * *c1;
	    *c0 = *c0 * 2.f * *c0 - 1.f;
	    *x *= 4.f;
	}
    }
    return 0;
} /* drift_kepu_stumpff__ */

/* ------------------------------------------------------------------ */

/* ************************************************************************* */
/*                        DRIFT_ONE.F */
/* ************************************************************************* */
/* This subroutine does the danby-type drift for one particle, using */
/* appropriate vbles and redoing a drift if the accuracy is too poor */
/* (as flagged by the integer iflg). */

/*             Input: */
/*                 nbod          ==>  number of massive bodies (int scalar) */
/*                 mass          ==>  mass of bodies (real array) */
/*                 x,y,z         ==>  initial position in jacobi coord */
/*                                    (real scalar) */
/*                 vx,vy,vz      ==>  initial position in jacobi coord */
/*                                    (real scalar) */
/*                 dt            ==>  time step */
/*             Output: */
/*                 x,y,z         ==>  final position in jacobi coord */
/*                                       (real scalars) */
/*                 vx,vy,vz      ==>  final position in jacobi coord */
/*                                       (real scalars) */
/*                 iflg          ==>  integer (zero for successful step) */

/* Authors:  Hal Levison & Martin Duncan */
/* Date:    2/10/93 */
/* Last revision: 2/10/93 */

/*   drift_kepu_stumpff */
/* Subroutine */ int drift_one__(doublereal *mu, doublereal *x, doublereal *y,
	 doublereal *z__, doublereal *vx, doublereal *vy, doublereal *vz, 
	doublereal *dt, integer *iflg)
{
    static integer i__;
    extern /* Subroutine */ int drift_dan__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal dttmp;

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
/* ...  Inputs and Outputs: */
/* ...  Output */
/* ...  Internals: */
/* ---- */
/* ...  Executable code */
    drift_dan__(mu, x, y, z__, vx, vy, vz, dt, iflg);
    if (*iflg != 0) {
	for (i__ = 1; i__ <= 10; ++i__) {
	    dttmp = *dt / 10.;
	    drift_dan__(mu, x, y, z__, vx, vy, vz, &dttmp, iflg);
	    if (*iflg != 0) {
		return 0;
	    }
	}
    }
    return 0;
} /* drift_one__ */

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
/* drift_one */
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
    static cilist io___1138 = { 0, 6, 0, 0, 0 };


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
    s_wsle(&io___1138);
    do_lio(&c__9, &c__1, "FGET : RETURNING WITHOUT COMPLETE CONVERGENCE", (
	    ftnlen)45);
    e_wsle();
    return ret_val;
} /* orbel_fget__ */

/* ------------------------------------------------------------------ */

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
/* orbel_fget */
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
/* orbel_fhybrid */
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
    static cilist io___1155 = { 0, 6, 0, 0, 0 };
    static cilist io___1157 = { 0, 6, 0, 0, 0 };
    static cilist io___1158 = { 0, 6, 0, 0, 0 };


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
    biga = pow_dd(&d__1, &c_b95);
    d__1 = b * .5f + sq;
    bigb = -pow_dd(&d__1, &c_b95);
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
    s_wsle(&io___1155);
    do_lio(&c__9, &c__1, "FLON : RETURNING WITHOUT COMPLETE CONVERGENCE", (
	    ftnlen)45);
    e_wsle();
    diff = *e * sinh(ret_val) - ret_val - *capn;
    s_wsle(&io___1157);
    do_lio(&c__9, &c__1, "N, F, ecc*sinh(F) - F - N : ", (ftnlen)28);
    e_wsle();
    s_wsle(&io___1158);
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
	tmp = pow_dd(&x, &c_b1082);
	ret_val = tmp - 1. / tmp;
    }
    if (iflag == 1) {
	ret_val = -ret_val;
	*q = -(*q);
    }
    return ret_val;
} /* orbel_zget__ */


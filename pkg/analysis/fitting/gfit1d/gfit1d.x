include	<imhdr.h>
include <pkg/gtools.h>
include	<error.h>
include	<math/curfit.h>
include	"icfit.h"
include	<mach.h>
include	<ctype.h>
include	<fset.h>	# FIO
include	<gset.h>	# GIO
include <tbset.h>	# TBTABLES

include "../lib/colnames.h"
include "../lib/curfitdef.h"
include "../lib/dcurfitdef.h"
include "../nlfit/nlfit.h"


define	SZ_BUF		2048		# Initial pixel buffer size
define	LIST_OP		1
define	IMAGE_OP	2
define	TABLE_OP	3
define	RE_TRY		-99		# error flag value
define	DEF_ERR		1.0		# default error bar size

# GFIT1D -- 1D fitting utility where input may be one or more lists (y or 
# x,y), image sections, or table columns.  Multidimensional image
# sections are reduced to a vector by computing the projection about the
# indicated axis.  The fitting parameters may be set interactively using the 
# icfit package
#
#
# 18-Jun-1990	I.Busko		Modified output from ASCII database 
#				to SDAS table, including coefficient errors
# 25-Jun-1990	I.Busko		Relinked with modified rdcurves routines, to
#				get info from image headers. Output of
#				power-series polynomial is supported.
# 13-Aug-1990	I.Busko		Double precision version
# 19-Sep-1997   I.Busko		Reset vector pointers to NULL
# 14-Sep-1998   I.Busko		Test for dummy 1st row in output table

procedure t_gfit1d ()

char	input[SZ_LINE]			# Input data names
char	output[SZ_FNAME]		# Output fit file
bool	interactive			# Interactive?
char	graphics[SZ_FNAME]		# Graphics device
char	sample[SZ_LINE]			# Sample ranges
int	naverage			# Sample averaging size
char	function[SZ_LINE]		# Curve fitting function
int	order				# Order of curve fitting function
real	low_reject, high_reject		# Rejection thresholds
int	niterate			# Number of rejection iterations
real	grow				# Rejection growing radius
int	axis				# Axis for projection
bool	ps				# Power-series output ?
char	errtyp[SZ_LINE]			# Source of error info
real	sigma				# Constant error bar
real	epadu, rerr			# ccd noise model parameters
real	xmin, xmax			# Normalization interval

pointer	xpt				# Pointer to x vector
pointer	ypt				# Pointer to y vector
pointer	size				# Pointer to error vector
pointer	out				# output table
pointer	ic				# ICFIT pointer
pointer	gt				# GTOOLS pointer
pointer	nl				# NLFIT pointer
pointer	ft, fty
char	word[SZ_FNAME,4]		# Input line words
char	image[SZ_FNAME]			# Image name
int	npix				# Number of values per vector
int	istat				# errors getting data
int	etype
long	cpu, clock			# time variables

int	gt_init(), strdic()
int	nnames, rc_rdcurves(), fstati()
int	clgeti()
real	clgetr()
long	cputime(), clktime()
bool	clgetb()
pointer	opendb()

begin
	# Get input and output files.

	if (fstati (STDIN, F_REDIR) == YES)
	    call strcpy ("STDIN", input, SZ_FNAME)
	else
	    call clgstr ("input", input, SZ_LINE)

	call clgstr ("output", output, SZ_FNAME)
	out = opendb (output)

	# Get task parameters.

	call clgstr ("sample", sample, SZ_LINE)
	naverage 	= clgeti ("naverage")
	call clgstr ("function", function, SZ_LINE)
	order 		= clgeti ("order")

	low_reject 	= clgetr ("low_reject")
	high_reject 	= clgetr ("high_reject")
	niterate 	= clgeti ("niterate")
	grow 		= clgetr ("grow")
	interactive 	= clgetb ("interactive")
	axis 		= clgeti ("axis")
	ps	 	= clgetb ("ps")
	xmin		= clgetr ("xmin")
	xmax		= clgetr ("xmax")
	sigma		= clgetr ("sigma")
	epadu		= clgetr ("epadu")
	rerr		= clgetr ("readnoise")
	call clgstr ("errtyp", errtyp, SZ_LINE)
	etype = strdic (errtyp, errtyp, SZ_LINE, NLERRORS)
	if (etype == 0)
	    call error (0, "Incorrect errtype specification.")

	if ( interactive )
	    call clgstr ("device", graphics, SZ_FNAME)

 	# initialize the curve fitting package.

	# Set the ICFIT pointer structure.
	call ic_open (ic)
	call ic_pstr (ic, "sample", sample)
	call ic_puti (ic, "naverage", naverage)
	call ic_pstr (ic, "function", function)
	call ic_puti (ic, "order", order)
	call ic_putr (ic, "low", low_reject)
	call ic_putr (ic, "high", high_reject)
	call ic_puti (ic, "niterate", niterate)
	call ic_putr (ic, "grow", grow)
	IC_FITERROR(ic) = NO
	gt = gt_init()

	# Initialize nlfit package. This is needed only in order
	# to pass error info to underlying prep_data routine.
	nl = NULL
	call nl_init (nl, 1, 0., true, 1, 2)
	call nl_putr (nl, "sigma", sigma)
	call nl_putr (nl, "epadu", epadu)
	call nl_putr (nl, "readnoise", rerr)
	call nl_puti (nl, "errtyp", etype)

	# change the default from "line" to "mark".  JC Hsu 6/15/94
	call gt_sets (gt, GTTYPE, "mark")

	# read the data one at a time
	ft   = NULL
	fty  = NULL
	xpt  = NULL
	ypt  = NULL
	size = NULL
	while ( rc_rdcurves( input, image, ft, fty, nnames, word, xpt, ypt, 
			     size, npix, axis, true, istat) != EOF ) {

	    if ( istat == RE_TRY)
	  	next

	    call ic_pstr (ic, "ylabel", image)

	    # Do it

	    cpu   = cputime (0)
	    clock = clktime (0)

	    iferr( call f1_gfit1d ( Memr[xpt], Memr[ypt], Memr[size], npix, 
		          	    ic, gt, nl, image, axis, interactive, 
	                            graphics, out, ps, xmin, xmax)) {
		call erract( EA_WARN)
		next
	    }

	    call printf ("%7.2f  CPU seconds,  %7.2f  elapsed minutes.\n")
	        call pargr (real (cputime (cpu)) / 1000.)
	        call pargr (real (clktime (clock)) / 60.)

	}

	call nl_free (nl)
	call ic_closed (ic)
	call gt_free (gt)
	call tbtclo( out)

end

# F1_GFIT1D -- Given the x, y and error data determine the fitting function.
# If the interactive flag is set then set the fitting parameters interactively.

procedure f1_gfit1d (x, y, err, npix, ic, gt, nl, filename, axis, interactive, 
		     graphics, out, ps, xmin, xmax)

real	x[ARB]				# input x data
real	y[ARB]				# input y data
real	err[ARB]			# input err data
int	npix				# number of pixels
pointer	ic				# ICFIT pointer
pointer	gt				# GTOOLS pointer
pointer	nl				# NLFIT pointer
char	filename[ARB]			# filename
int	axis				# image axis
bool	interactive			# Interactive?
char	graphics[SZ_FNAME]		# Graphics device
pointer	out				# database pointer
bool	ps				# power-series output >
real	xmin, xmax			# normalization interval

int compare() # function to be called to compare elements
extern	compare

int	i
real	xmin1, xmax1
pointer	cv, gp, sp
pointer	ttx, tx, ty, twts
pointer	dx, dy, dw, wts, indpt, gopen()

common	/comm/ttx

errchk	gopen()

begin
	# Allocate memory for curve fitting.
	call smark (sp)
	call salloc (wts, npix, TY_REAL)
	call salloc (tx , npix, TY_REAL)
	ttx	= tx
	call salloc (ty , npix, TY_REAL)
	call salloc (twts , npix, TY_REAL)
	call salloc (indpt, npix, TY_INT)

	# Prepare weight vector.
	call prep_data (x, y, err, npix, nl, wts)

	# Sort the data for use before fitting.
	do i = 1, npix
	   Memi[indpt+i-1] = i
	call amovr( x, Memr[tx], npix)
	call tqsort ( Memi[indpt], npix, compare)
	do i = 1, npix {
	    Memr[tx+i-1] = x[ Memi[indpt+i-1] ]
	    Memr[ty+i-1] = y[ Memi[indpt+i-1] ]
	    Memr[twts+i-1] = Memr[wts+ Memi[indpt+i-1] - 1]
	}

	# Find normalization interval.
	xmin1 = xmin
	xmax1 = xmax
	if (IS_INDEFR (xmin))
	    xmin1 = Memr[tx]
	if (IS_INDEFR (xmax))
	    xmax1 = Memr[tx+npix-1]
	call ic_putr (ic, "xmin", xmin1 )
	call ic_putr (ic, "xmax", xmax1 )

	# Transfer data to double precision vectors.
	call salloc (dx, npix, TY_DOUBLE)
	call salloc (dy, npix, TY_DOUBLE)
	call salloc (dw, npix, TY_DOUBLE)
	do i = 0, npix-1 {
	    Memd[dx+i] = Memr[tx+i]
	    Memd[dy+i] = Memr[ty+i]
	    Memd[dw+i] = Memr[twts+i]
	}

	# If the interactive flag is set then use icg_fit to set the
	# fitting parameters.  
	# The weights are reset since the user may delete points.

	if (interactive) {
	    gp = gopen (graphics, NEW_FILE, STDGRAPH)

	    call icg_fitd (ic, gp, "cursor", gt, cv, 
				Memd[dx], Memd[dy], Memd[dw], npix)
	    call gclose (gp)
	} else {
	    iferr (call ic_fitd (ic, cv, Memd[dx], Memd[dy], Memd[dw], npix, 
				 YES, YES, YES, YES ))
		IC_FITERROR(ic) = YES
	}

	call f1d_tdb( filename, ic, cv, out, ps, Memd[dx], Memd[dy], Memd[dw],
	              npix)

	call dcvfree (cv)
	call sfree (sp)

end

# F1D_TDB -- write fit details to table

procedure f1d_tdb( filename, ic, cv, tp, ps, x, y, wts, npts )

char	filename[ARB]	# i: filename
pointer	ic		# i: icfit pointer
pointer	cv		# i: curfit pointer
pointer	tp		# i: table pointer
bool	ps		# i: power-series output ?
double	x[ARB]		# i: independent variable
double	y[ARB]		# i: dependent variable
double	wts[ARB]	# i: weight
int	npts		# i: total number of data points

char	str[SZ_LINE], func[SZ_LINE]
int	row, tnpar, buflen, ncoeff, ifunc, i, n, deleted
long	time
real	chisq, rms
double	aux
pointer	sp, coeff, colptr, fit, error, wts1
bool	null

int	dcvstati(), tbpsta()
long	clktime()
real	ic_getr()
bool	streq()

begin
	call smark( sp)
	buflen	= dcvstati( cv, CVNSAVE)
	call salloc( coeff, buflen, TY_DOUBLE )

	# Number of coefficients
	ncoeff = buflen - 4

	# If necessary, increase number of columns in table
	tnpar = (tbpsta (tp, TBL_NCOLS) - MIN_COLS) / 2
	if (tnpar < ncoeff) {
	    do i = tnpar + 1, ncoeff {
	        call sprintf (str, SZ_LINE, DB_CCOEF)
	            call pargi (i)
	        call tbcdef (tp, colptr, str, DB_UCOEF, DB_FCOEF, TY_REAL, 1, 1)
	        call sprintf (str, SZ_LINE, DB_CERR)
	            call pargi (i)
	        call tbcdef (tp, colptr, str, DB_UERR, DB_FERR, TY_REAL, 1, 1)
	    }
	}

	# This code tests the contents of the first row in the table,
        # to see if it is a dummy row created by routine opendb, or has
        # valid contents. A dummy row is overwritten. This was put in
        # place to correct the behavior of the FITSIO library. It crashesed
        # the program when attempting to add new columns to a table that
        # has no rows (9/14/98, IB)
	row = tbpsta (tp, TBL_NROWS)
	if (row == 1) {
	    call tbcfnd (tp, DB_CFILE, colptr, 1)
	    call tbrgtt (tp, colptr, str, null, DB_SFILE, 1, 1)
	    if (!streq (str, DUMMY))
	        row = 2
	} else
	    # If more than one row in table, just append.
	    row = row + 1

	call tbcfnd (tp, DB_CFILE, colptr, 1)	# write file name
	call tbeptt (tp, colptr, row, filename)

	time = clktime (0)			# write time stamp
	call cnvtime (time, str, SZ_LINE)
	call tbcfnd (tp, DB_CTIME, colptr, 1)
	call tbeptt (tp, colptr, row, str)

	# get coefficients from curfit
	call dcvsave( cv, Memd[coeff] )

	# set function
	ifunc = int(Memd[coeff])
	if (ps && ((ifunc == CV_CHEB) || (ifunc == CV_LEG)))
	    call strcpy (PSPOLY, func, SZ_LINE)
	else
	    call extnstr (FUNCTIONS, ifunc, func) 
	call sprintf (str, SZ_LINE, "%s")
	    call pargstr (func)
	call tbcfnd (tp, DB_CFUNC, colptr, 1)
	call tbeptt (tp, colptr, row, str)

	call tbcfnd (tp, DB_CUNIT, colptr, 1)		# unit
	call tbeptt (tp, colptr, row, "*")

	call tbcfnd (tp, DB_CDEGR, colptr, 1)		# no. of coefficients
	call tbepti (tp, colptr, row, ncoeff)

	call tbcfnd (tp, DB_CNPTS, colptr, 1)		# npts
	call tbepti (tp, colptr, row, npts)

	call tbcfnd (tp, DB_CXMIN, colptr, 1)		# xmin, xmax
	call tbeptr (tp, colptr, row, ic_getr( ic, "xmin"))
	call tbcfnd (tp, DB_CXMAX, colptr, 1)
	call tbeptr (tp, colptr, row, ic_getr( ic, "xmax"))

	# obtain power-series equivalents
	if (ps) {
	    if (streq (func, PSPOLY))
#	        call dcvpower (cv, Memd[coeff+4], buflen-4)
# This is a temporary fix until IRAF VMS is upgraded to 10.4
	        call t_dcvpower (cv, Memd[coeff+4], buflen-4)
	}

	do i = 5, buflen {				# coefficients
	    call sprintf (str, SZ_LINE, DB_CCOEF)
	        call pargi (i-4)
	    call tbcfnd (tp, str, colptr, 1)
	    call tbeptd (tp, colptr, row, Memd[coeff+i-1])
	}

	# get errors directly from icfit
	call salloc( error, ncoeff, TY_DOUBLE )
	n = IC_NFIT(ic)
	deleted = 0
	chisq = INDEF
	rms   = INDEF
	if (n == npts) {
	    call salloc (fit, n, TY_DOUBLE)
	    call salloc (wts1, n, TY_DOUBLE)
	    # eliminate rejected points and count deleted points.
	    call amovd (wts, Memd[wts1], n)
	    if (IC_NREJECT(ic) > 0) {
		do i = 1, npts {
		    if (Memi[IC_REJPTS(ic)+i-1] == YES)
			Memd[wts1+i-1] = 0.
		}
	    }
	    deleted = 0
	    do i = 1, n {
		if (wts[i] == 0.)
		    deleted = deleted + 1
	    }
	    # set the fit and compute the RMS error.
	    if (IC_FITERROR(ic) == NO) {
	        call dcvvector (cv, x, Memd[fit], n)
	        if (ps && streq (func, PSPOLY))
	            call dcvepower (cv, y, wts, Memd[fit], npts, aux, 
	                            Memd[error])
	        else
	            call dcverrors (cv, y, wts, Memd[fit], npts, aux, 
	                            Memd[error])
	        call f1d_rms (x, y, Memd[fit], Memd[wts1], npts, ncoeff,
	                      chisq, rms)
	    } else {
		rms   = INDEF
	        chisq = INDEF
	    }
	} else if (n > 0) {
	    # Allocate memory for the fit.
	    call salloc (fit, n, TY_DOUBLE)
	    call salloc (wts1, n, TY_DOUBLE)
	    # eliminate rejected points and count deleted points.
	    call amovd (Memd[IC_WTSFIT(ic)], Memd[wts1], n)
	    if (IC_NREJECT(ic) > 0) {
	        do i = 1, npts {
		    if (Memi[IC_REJPTS(ic)+i-1] == YES)
			Memd[wts1+i-1] = 0.
		}
	    }
	    deleted = 0
	    do i = 1, n {
		if (wts[i] == 0.)
		    deleted = deleted + 1
	    }
	    # set the fit and compute chisq and rms.
	    if (IC_FITERROR(ic) == NO) {
	        call dcvvector (cv, Memd[IC_XFIT(ic)], Memd[fit], n)
	        if (ps && streq (func, PSPOLY))
	            call dcvepower (cv, Memd[IC_YFIT(ic)], Memd[wts1], 
                                    Memd[fit], n, aux, Memd[error])
	        else
	            call dcverrors (cv, Memd[IC_YFIT(ic)], Memd[wts1], 
                                    Memd[fit], n, aux, Memd[error])
	        call f1d_rms (Memd[IC_XFIT(ic)], Memd[IC_YFIT(ic)], Memd[fit], 
	                      Memd[wts1], n, ncoeff, chisq, rms)
	    } else {
	        chisq = INDEF
		rms   = INDEF
	    }
	}

	call tbcfnd (tp, DB_CRMS, colptr, 1)		# rms
	call tbeptr (tp, colptr, row, rms)

	call tbcfnd (tp, DB_CCHI, colptr, 1)		# chi-sq
	call tbeptr (tp, colptr, row, chisq)

	# Add 6/27/94 by JC Hsu, to properly calculate the errors of the
	# power series coefficients
#	if (ps) {
#	    if (streq (func, PSPOLY))
#	        call cvpserr (cv, Memd[error], buflen-4)
#	}
# The above was commented out and replaced by the calls to curfit's dcvepower
# just for the sake of consistency; both routines generate the same results.
# Detailed analysis of error computation will be published elsewhere 
# (IB 06/03/96).

	do i = 1, ncoeff {				# errors
	    call sprintf (str, SZ_LINE, DB_CERR)
	        call pargi (i)
	    call tbcfnd (tp, str, colptr, 1)
	    call tbeptd (tp, colptr, row, Memd[error+i-1])
	}

	call sfree (sp)

	call printf ("Line %d created in output table.\n")
	call pargi (row)
end


# F1D_RMS -- Compute chisq and rms of points which have not been deleted.

procedure f1d_rms (x, y, fit, wts, npts, ncoeff, chisq, rms)

double	x[ARB]		# i: Independent variable
double	y[ARB]		# i: Dependent variable
double	fit[ARB]	# i: Fit
double	wts[ARB]	# i: Weights. Assume they are the inverse of the
			#    y errors.
int	npts		# i: Number of data points
int	ncoeff		# i: Number of coefficients
real	chisq		# o: Chi-square
real	rms		# o: rms

int	i, n
double	resid, schi, srms

begin
	schi = 0.d0
	srms = 0.d0
	n = 0
	do i = 1, npts {
	    if (wts[i] == 0.)
		next
	    resid = y[i] - fit[i]
	    srms  = srms + resid * resid
	    resid = resid * wts[i]
	    schi  = schi + resid * resid
	    n = n + 1
	}

	if (n > 0)
	    rms = real(sqrt (srms / n))
	else
	    rms = INDEF

	if (n > ncoeff)
	    chisq = real (schi / (n - ncoeff))
	else
	    chisq = INDEF
end

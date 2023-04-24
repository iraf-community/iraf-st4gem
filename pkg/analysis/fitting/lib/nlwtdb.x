include	<time.h>
include	<error.h>
include	<mach.h>
include	<ctype.h>
include <fset.h>
include <tbset.h>

include "../nlfit/nlfit.h"
include	"colnames.h"

# NL_WTDB -- write non-linear fit details to database table
#
# Modified 2/22/96 (I Busko) to implement chi-squared minimization.
# Modified 9/14/98 (IB) to write dummy row at table creation time. This
#    enables later creation of new columns in FITS tables.

procedure nl_wtdb (filename, nl, tp, verbose)

char	filename[ARB]	# original file name
pointer	nl		# nlfit pointer.
pointer	tp		# table pointer.
bool	verbose		# print results ?

char	str[SZ_LINE], str1[SZ_LINE]
int	npts, row, tnpar, buflen, funct, unit, i
long	time
real	rms, chisq
pointer	sp, coeff, error, colptr
bool	null

int	nl_stati(), tbpsta()
long	clktime()
real	nl_statr()
bool	streq()

begin
	# get fit parameters
	call smark( sp)
	buflen = nl_stati( nl, "npar")
	call salloc( coeff, buflen, TY_REAL )
	call salloc( error, buflen, TY_REAL )
	call nl_gcoeff( nl, Memr[coeff], buflen )	# coefficients
	call nl_gerrors( nl, Memr[error], buflen )	# errors
	rms   = nl_statr (nl, "rms")
	chisq = nl_statr (nl, "chisq")
	funct = nl_stati (nl, "fitfunc")
	unit  = nl_stati (nl, "units")
	npts  = nl_stati (nl, "npts") - nl_stati (nl, "nreject")

	if (tp != NULL) {
	    # if necessary, increase number of columns in table
	    tnpar = (tbpsta (tp, TBL_NCOLS) - MIN_COLS) / 2
	    if (tnpar < buflen) {
	        do i = tnpar + 1, buflen {
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

	    call tbcfnd (tp, DB_CFILE, colptr, 1)	# file name
	    call tbeptt (tp, colptr, row, filename)

	    time = clktime (0)			# time
	    call cnvtime (time, str, SZ_LINE)
	    call tbcfnd (tp, DB_CTIME, colptr, 1)
	    call tbeptt (tp, colptr, row, str)

	    call extnstr (NLFUNCTIONS, funct, str1)	# function
	    call sprintf (str, SZ_LINE, "%s")
	        call pargstr (str1)
	    call tbcfnd (tp, DB_CFUNC, colptr, 1)
	    call tbeptt (tp, colptr, row, str)

	    call tbcfnd (tp, DB_CDEGR, colptr, 1)	# no. of coefficients
	    call tbepti (tp, colptr, row, buflen)

	    call extnstr (NLUNITS, unit, str1)	# unit
	    call sprintf (str, SZ_LINE, "%s")
	        call pargstr (str1)
	    call tbcfnd (tp, DB_CUNIT, colptr, 1)
	    call tbeptt (tp, colptr, row, str)

	    call tbcfnd (tp, DB_CNPTS, colptr, 1)	# npts
	    call tbepti (tp, colptr, row, npts)

	    call tbcfnd (tp, DB_CRMS, colptr, 1)	# rms
	    call tbeptr (tp, colptr, row, rms)

	    call tbcfnd (tp, DB_CCHI, colptr, 1)	# chisq
	    call tbeptr (tp, colptr, row, chisq)

	    do i = 1, buflen {
	        call sprintf (str, SZ_LINE, DB_CCOEF)
	            call pargi (i)
	        call tbcfnd (tp, str, colptr, 1)
	        call tbeptr (tp, colptr, row, Memr[coeff+i-1])
	    }

	    do i = 1, buflen {
	        call sprintf (str, SZ_LINE, DB_CERR)
	            call pargi (i)
	        call tbcfnd (tp, str, colptr, 1)
	        call tbeptr (tp, colptr, row, Memr[error+i-1])
	    } 
	}

	if (verbose) {					# print at STDOUT
	    switch (funct) {
	    case TWODGAUSS:
	        call printf ("\nBackg = %g  (%g)")
	            call pargr (Memr[coeff+NL_G2A])
	            call pargr (Memr[error+NL_G2A])
	        call printf ("\tChi-sq = %g\n")
	            call pargr (chisq)
	        call printf ("\tResid = %g\n")
	            call pargr (rms)
	        call printf ("Ampl  = %g  (%g)\n")
	            call pargr (Memr[coeff+NL_G2AMPL])
	            call pargr (Memr[error+NL_G2AMPL])
	        call printf ("Xcent = %g  (%g)")
	            call pargr (Memr[coeff+NL_G2XC])
	            call pargr (Memr[error+NL_G2XC])
	        call printf ("\tYcent = %g  (%g)\n")
	            call pargr (Memr[coeff+NL_G2YC])
	            call pargr (Memr[error+NL_G2YC])
	        call printf ("FWHM  = %g  (%g)\n")
	            call pargr (Memr[coeff+NL_G2FWHM])
	            call pargr (Memr[error+NL_G2FWHM])
	        call printf ("Ellip = %g  (%g)")
	            call pargr (Memr[coeff+NL_G2ELL])
	            call pargr (Memr[error+NL_G2ELL])
	        call printf ("\tTheta  = %g  (%g)\n")
	            call pargr (Memr[coeff+NL_G2TETA]/3.1415926*180.)
	            call pargr (Memr[error+NL_G2TETA]/3.1415926*180.)


	    default:

	    }

	    if (tp != NULL) {
	        call printf ("Line %d created in output table.\n")
	            call pargi (row)
	    }
	}
	call sfree (sp)
end

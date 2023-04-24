include <tbset.h>

include "colnames.h"

# GETFIT  --  Reads one line from a fit table. Arrays for the coefficients, 
# errors and flags are allocated on the stack. Caller must set stack marker.
#
# Modified 2/22/96 (I Busko) to implement chi-squared handling.

procedure getfit (table, row, origfile, timestamp, function, unit, strsize,
		  npar, npts, rms, chisq, xmin, xmax, coeff, err, coefflags)

char	table[ARB]			# i: fit table
int	row				# i: table row
char	origfile[strsize]		# o: original file name
char	timestamp[strsize]		# o: time stamp
char	function[strsize]		# o: function type
char	unit[strsize]			# o: units
int	strsize				# i: string size
int	npar				# o: number of coefficients 
int	npts				# o: number of data points
real	rms				# o: rms of fit
real	chisq				# o: chisq of fit
real	xmin, xmax			# o: extremes for normalization
pointer coeff				# o: coefficient vector
pointer	err				# o: errors vector
pointer coefflags			# o: flags vector

char	str[SZ_LINE]
int	i
pointer	dt, colptr
bool	null

int	tbpsta()
pointer	opendb()

begin
	# open table
	dt = opendb (table)
	if (row > tbpsta (dt, TBL_NROWS))
	    call error (0, "Error in row specification.")

	# get original file name
	call tbcfnd (dt, DB_CFILE, colptr, 1)
	call tbrgtt (dt, colptr, origfile, null, strsize, 1, row)

	# get time stamp
	call tbcfnd (dt, DB_CTIME, colptr, 1)
	call tbrgtt (dt, colptr, timestamp, null, strsize, 1, row)

	# get function type
	call tbcfnd (dt, DB_CFUNC, colptr, 1)
	call tbrgtt (dt, colptr, function, null, strsize, 1, row)

	# get units
	call tbcfnd (dt, DB_CUNIT, colptr, 1)
	call tbrgtt (dt, colptr, unit, null, strsize, 1, row)

	# get npts
	call tbcfnd (dt, DB_CNPTS, colptr, 1)
	call tbrgti (dt, colptr, npts, null, 1, row)

	# get rms
	call tbcfnd (dt, DB_CRMS, colptr, 1)
	call tbrgtr (dt, colptr, rms, null, 1, row)

	# get chi-sq
	call tbcfnd (dt, DB_CCHI, colptr, 1)
	call tbrgtr (dt, colptr, chisq, null, 1, row)

	# get xmin, xmax
	call tbcfnd (dt, DB_CXMIN, colptr, 1)
	call tbrgtr (dt, colptr, xmin, null, 1, row)
	call tbcfnd (dt, DB_CXMAX, colptr, 1)
	call tbrgtr (dt, colptr, xmax, null, 1, row)

	# alloc space for coefficients
	call tbcfnd (dt, DB_CDEGR, colptr, 1)
	call tbrgti (dt, colptr, npar, null, 1, row)
	call salloc ( coeff, npar, TY_REAL )
	call salloc ( err, npar, TY_REAL)
	call salloc ( coefflags, npar, TY_BOOL)

	# read coefficients
	do i = 1, npar {
	    call sprintf( str, SZ_LINE, DB_CCOEF)
		call pargi( i)
	    call tbcfnd (dt, str, colptr, 1)
	    call tbrgtr (dt, colptr, Memr[coeff+i-1], null, 1, row)
	}

	# read errors
	do i = 1, npar {
	    call sprintf( str, SZ_LINE, DB_CERR)
	        call pargi( i)
	    call tbcfnd (dt, str, colptr, 1)
	    call tbrgtr (dt, colptr, Memr[err+i-1], null, 1, row)
	}

	# set flag array
	# do i = 0, npar-1
	#   Memb[coefflags+i] = false

	# set the flags according to the error values, if less than 0,
	# set to false, i.e. does not vary. Taken from ngaussfit.
	do i = 1, npar {
	    Memb[coefflags+i-1] = (Memr[err+i-1] >= 0.)
	}

	call tbtclo (dt)
end

include <math/gsurfit.h>

# xyz_fit -- fit the function to the data
# The gsurfit package is used to fit a function to the (x, y, z) data,
# and the sf pointer is returned.  The calling routine should release
# sf when no longer needed.
#
# Phil Hodge, 15-Dec-1993  Subroutine created.

procedure xyz_fit (sf, x, y, z, nrows, ndim,
		func, xorder, yorder, cross_terms, x1, x2, y1, y2, rms)

pointer sf		# o: pointer to gsurfit struct
double	x[ARB]		# i: X values
double	y[ARB]		# i: Y values (if 2-D)
double	z[ARB]		# i: Z values
int	nrows		# i: size of x, y, z arrays
int	ndim		# i: dimension of fit, 1 or 2
int	func		# i: indicates which function to fit
int	xorder, yorder	# i: order (plus one) of fit in X and Y
bool	cross_terms	# i: true if cross terms are to be included
double	x1, x2		# i: limits on X values
double	y1, y2		# i: limits on Y values
double	rms		# o: RMS deviation of fit from input
#--
double	sum, sumsq, xn	# for getting RMS
double	v		# value of fitted function at a point
double	weight		# returned by gsaccum and then ignored
int	xterms		# include cross terms (YES or NO)?
int	ier		# error return code from gssolve
int	i
double	dgseval()

begin
	if (cross_terms)
	    xterms = YES
	else
	    xterms = NO

	# Initialize gsurfit.
	call dgsinit (sf, func, xorder, yorder, xterms, x1, x2, y1, y2)

	# Accumulate sums.
	if (ndim == 1) {
	    do i = 1, nrows
		call dgsaccum (sf, x[i], 1.d0, z[i], weight, WTS_UNIFORM)
	} else {			# 2-D
	    do i = 1, nrows
		call dgsaccum (sf, x[i], y[i], z[i], weight, WTS_UNIFORM)
	}

	# Do the fit.
	call dgssolve (sf, ier)

	if (ier == SINGULAR)
	    call error (1, "matrix is singular")
	else if (ier == NO_DEG_FREEDOM)
	    call error (1, "not enough data for specified order of fit")
	else if (ier != OK)
	    call error (1, "unspecified error return from dgssolve")

	if (nrows < 2) {			# can't compute RMS
	    rms = 0.d0
	    return
	}

	# Find the RMS deviation of the fit from the input values.
	sum = 0.d0
	sumsq = 0.d0
	if (ndim == 1) {			# 1-D
	    do i = 1, nrows {
		# Evaluate the fit at x[i], and subtract actual value.
		v = dgseval (sf, x[i], 1.d0) - z[i]
		sum = sum + v
		sumsq = sumsq + v*v
	    }
	} else {				# 2-D
	    do i = 1, nrows {
		# Evaluate the fit at x[i],y[i], and subtract actual value.
		v = dgseval (sf, x[i], y[i]) - z[i]
		sum = sum + v
		sumsq = sumsq + v*v
	    }
	}

	xn = double(nrows)
	rms = sqrt ((sumsq - sum*sum/xn) / (xn - 1.d0))
end

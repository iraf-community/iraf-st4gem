# Copyright(c) 1986 Association of Universities for Research in Astronomy Inc.

include <mach.h>
include	<math/curfit.h>

include	"dcurfitdef.h"

# CVPOWER -- Convert legendre or chebyshev coeffecients to power series.

procedure t_dcvpower (cv, ps_coeff, ncoeff)

pointer	cv				# Pointer to curfit structure
double	ps_coeff[ncoeff]		# Power series coefficients (output)
int	ncoeff				# Number of coefficients in fit

pointer	sp, cf_coeff, elm
int	function
int	dcvstati()

begin
	function = dcvstati (cv, CVTYPE)

	if (function != LEGENDRE && function != CHEBYSHEV) {
	    call eprintf ("Cannot convert coefficients - wrong function type\n")
	    call amovkd (INDEFD, ps_coeff, ncoeff)
	    return
	}

	call smark (sp)
	call salloc (elm, ncoeff ** 2, TY_DOUBLE)
	call salloc (cf_coeff, ncoeff, TY_DOUBLE)

	call amovkd (0.0d0, Memd[elm], ncoeff ** 2)

	# Get existing coefficients
	call dcvcoeff (cv, Memd[cf_coeff], ncoeff)

	switch (function){
	case (LEGENDRE):
	    call dcv_mlegen (Memd[elm], ncoeff)
	    call dcv_legen (Memd[elm], Memd[cf_coeff], ps_coeff, ncoeff)
	case (CHEBYSHEV):
	    call dcv_mcheby (Memd[elm], ncoeff)
	    call dcv_cheby (Memd[elm], Memd[cf_coeff], ps_coeff, ncoeff)
	}

	# Normalize coefficients
	call dcv_normalize (cv, ps_coeff, ncoeff)

	call sfree (sp)
end


# CVEPOWER -- Procedure to calculate the reduced chi-squared of the fit
# and the standard deviations of the power series coefficients. First the
# variance and the reduced chi-squared of the fit are estimated. If these
# two quantities are identical the variance is used to scale the errors
# in the coefficients. The errors in the coefficients are proportional
# to the inverse diagonal elements of MATRIX.

procedure dcvepower (cv, y, w, yfit, npts, chisqr, perrors)

pointer	cv		# curve descriptor
double	y[ARB]		# data points
double	yfit[ARB]	# fitted data points
double	w[ARB]		# array of weights
int	npts		# number of points
double	chisqr		# reduced chi-squared of fit
double	perrors[ARB]	# errors in coefficients

int	i, j, n, nfree, function, ncoeff
double	variance, chisq, hold
pointer	sp, covar, elm
int	dcvstati()

begin
	# Determine the function type.
	function = dcvstati (cv, CVTYPE)
	ncoeff = dcvstati (cv, CVNCOEFF)

	# Check the function type.
	if (function != LEGENDRE && function != CHEBYSHEV) {
	    call eprintf ("Cannot convert errors - wrong function type\n")
	    call amovkd (INDEFD, perrors, ncoeff)
	    return
	}

        # Estimate the variance and chi-squared of the fit.
        n = 0
        variance = 0.
        chisq = 0.
        do i = 1, npts {
            if (w[i] <= 0.0)
                next
            hold = (y[i] - yfit[i]) ** 2
            variance = variance + hold
            chisq = chisq + hold * w[i]
            n = n + 1
        }

        # Calculate the reduced chi-squared.
        nfree = n - CV_NCOEFF(cv)
        if (nfree  > 0)
            chisqr = chisq / nfree
        else
            chisqr = 0.

        # If the variance equals the reduced chi_squared as in the case of
	# uniform weights then scale the errors in the coefficients by the
	# variance not the reduced chi-squared
        if (abs (chisq - variance) <= DELTA) {
            if (nfree > 0)
                variance = chisq / nfree
            else
                variance = 0.
        } else
            variance = 1.


	# Allocate space for the covariance and conversion matrices.
	call smark (sp)
	call salloc (covar, ncoeff * ncoeff, TY_DOUBLE)
	call salloc (elm, ncoeff * ncoeff, TY_DOUBLE)

	# Compute the covariance matrix.
	do j = 1, ncoeff {
	    call aclrd (perrors, ncoeff)
	    perrors[j] = double(1.0)
	    call dcvchoslv (CHOFAC(CV_CHOFAC(cv)), CV_ORDER(cv),
	        CV_NCOEFF(cv), perrors, perrors)
	    call amulkd (perrors, double(variance), perrors, ncoeff)
	    call achtdd (perrors, Memd[covar+(j-1)*ncoeff], ncoeff)
	}

	# Compute the conversion matrix.
	call amovkd (0.0d0, Memd[elm], ncoeff * ncoeff)
	switch (function) {
	case LEGENDRE:
	    call dcv_mlegen (Memd[elm], ncoeff)
	case CHEBYSHEV:
	    call dcv_mcheby (Memd[elm], ncoeff)
	}

	# Normalize the errors to the appropriate data range.
	call dcv_enormalize (cv, Memd[elm], ncoeff)

	# Compute the new squared errors.
	call dcv_etransform (cv, Memd[covar], Memd[elm], perrors, ncoeff)

	# Compute the errors.
	do j = 1, ncoeff {
	    if (perrors[j] >= 0.0)
	        perrors[j] = sqrt(perrors[j])
	    else
		perrors[j] = 0.0
	}

	call sfree (sp)
end


# CV_MLEGEN -- Compute the matrix required to convert from legendre
# coefficients to power series coefficients. Summation notation for Legendre
# series taken from Arfken, page 536, equation 12.8.

procedure dcv_mlegen (matrix, ncoeff)

double	matrix[ncoeff, ncoeff]
int	ncoeff

int	s, n, r
double	dcv_legcoeff()

begin
	# Calculate matrix elements.
	do s = 0, ncoeff - 1 {
	    if (mod (s, 2) == 0) 
	        r = s / 2
	    else 
	        r = (s - 1) / 2

	    do n = 0, r
		matrix[s+1, (s+1) - (2*n)] = dcv_legcoeff (n, s)
	}
end


# CV_ETRANSFORM -- Convert the square of the fitted polynomial errors
# to the values appropriate for the equivalent power series polynomial.

procedure dcv_etransform (cv, covar, elm, perrors, ncoeff)

pointer	cv
double	covar[ncoeff,ncoeff]
double	elm[ncoeff,ncoeff]
double	perrors[ncoeff]
int	ncoeff

int	i, j, k
double	sum

begin
	do i = 1, ncoeff {
	    sum = 0.0d0
	    do j = 1, ncoeff {
		sum = sum + elm[j,i] * covar[j,j] * elm[j,i]
		do k = j + 1, ncoeff {
		    sum = sum + 2.0 * elm[j,i] * covar[j,k] * elm[k,i]
		}
	    }
	    perrors[i] = sum
	}
end


# CV_LEGEN -- Convert legendre coeffecients to power series coefficients.
# Scaling the coefficients from -1,+1 to the full data range is done in a 
# seperate procedure (cf_normalize).

procedure dcv_legen (matrix, cf_coeff, ps_coeff, ncoeff)

double	matrix[ncoeff, ncoeff]
double	cf_coeff[ncoeff]
double	ps_coeff[ncoeff]
int	ncoeff

int	n, i
double	sum

begin
	# Multiply matrix columns by curfit coefficients and sum.
	do n = 1, ncoeff {
	    sum = 0.0d0
	    do i = 1, ncoeff
	        sum = sum + (matrix[i,n] * cf_coeff[i])
	    ps_coeff[n] = sum
	}
end


# CV_LEGCOEFF -- calculate matrix elements for converting legendre coefficients
# to powers of x.

double procedure dcv_legcoeff (k, n)

int	k
int	n

double	fcn, sum1, divisor
double	dcv_factorial()

begin
	sum1 = ((-1) ** k) * dcv_factorial (2 * n - 2 * k)
	divisor = (2**n) * dcv_factorial (k) * dcv_factorial (n-k) * 
	    dcv_factorial (n - 2*k)
	fcn = sum1 / divisor

	return (fcn)
end


# CV_MCHEBY -- Compute the matrix required to convert from Chebyshev
# coefficient to power series coefficients. Summation notation for Chebyshev
# series from Arfken, page 628, equation 13.83

procedure dcv_mcheby (matrix, ncoeff)

double	matrix[ncoeff, ncoeff]		# Work array for matrix elements
int	ncoeff				# Number of coefficients

int	s, n, m
double	dcv_chebcoeff()

begin
	# Set first matrix element.
	matrix[1,1] = 1.0d0

	# Calculate remaining matrix elements.
	do s = 1, ncoeff - 1 {
	    if (mod (s, 2) == 0)
	        n = s / 2
	    else 
	        n = (s - 1) / 2

	    do m = 0, n
		matrix[(s+1),(s+1)-(2*m)] = (double(s)/2.0) *
		    dcv_chebcoeff (m, s)
	}
end


# CV_CHEBY -- Convert chebyshev coeffecients to power series coefficients.
# Scaling the coefficients from -1,+1 to the full data range is done in a 
# seperate procedure (cf_normalize).

procedure dcv_cheby (matrix, cf_coeff, ps_coeff, ncoeff)

double	matrix[ncoeff, ncoeff]		# Work array for matrix elements
double	cf_coeff[ncoeff]		# Input curfit coefficients
double	ps_coeff[ncoeff]		# Output power series coefficients
int	ncoeff				# Number of coefficients

int	n, i
double	sum

begin
	# Multiply matrix columns by curfit coefficients and sum.
	do n = 1, ncoeff {
	    sum = 0.0d0
	    do i = 1, ncoeff
	        sum = sum + (matrix[i,n] * cf_coeff[i])
	    ps_coeff[n] = sum
	}
end


# CV_CHEBCOEFF -- calculate matrix elements for converting chebyshev 
# coefficients to powers of x.

double procedure dcv_chebcoeff (m, n)

int	m	# Summation notation index
int	n	# Summation notation index

double	fcn, sum1, divisor
double	dcv_factorial()

begin
	sum1 = ((-1) ** m) * dcv_factorial (n - m - 1) * (2 ** (n - (2*m)))
	divisor = dcv_factorial (n - (2*m)) * dcv_factorial (m)
	fcn = sum1 / divisor

	return (fcn)
end


# CV_NORMALIZE -- Return coefficients scaled to full data range.

procedure dcv_normalize (cv, ps_coeff, ncoeff)

pointer	cv			# Pointer to curfit structure
int	ncoeff			# Number of coefficients in fit
double	ps_coeff[ncoeff]	# Power series coefficients

pointer	sp, elm, index
int	n, i, k
double	k1, k2, bc, sum

double	dcv_bcoeff()

begin
	# Need space for ncoeff**2 matrix elements
	call smark (sp)
	call salloc (elm, ncoeff ** 2, TY_DOUBLE)

	k1 = CV_RANGE(cv)
	k2 = k1 * CV_MAXMIN(cv)

	# Fill matrix, after zeroing it. 
	call amovkd (0.0d0, Memd[elm], ncoeff ** 2)
	do n = 1, ncoeff {
	    k = n - 1
	    do i = 0, k {
		bc = dcv_bcoeff (k, i)
		index = elm + k * ncoeff + i
		Memd[index] =  bc * ps_coeff[n] * (k1 ** i) * (k2 ** (k-i))
	    }
	}

	# Now sum along matrix columns to get coefficient of individual 
	# powers of x.
	do n = 1, ncoeff {
	   sum = 0.0d0
	   do i = 1, ncoeff {
	       index = elm + (n-1) + (i-1) * ncoeff
	       sum = sum + Memd[index]
	    }
	    ps_coeff[n] = sum
	}

	call sfree (sp)
end


# CV_ENORMALIZE -- Return the squares of the errors scaled to full data range.

procedure dcv_enormalize (cv, elm, ncoeff)

pointer	cv			# Pointer to curfit structure
double	elm[ncoeff,ncoeff]	# Input transformed matrix
int	ncoeff			# Number of coefficients in fit

pointer	sp, norm, onorm, index
int	n, i, k
double	k1, k2, bc

double	dcv_bcoeff()

begin
	# Need space for ncoeff**2 matrix elements
	call smark (sp)
	call salloc (norm, ncoeff ** 2, TY_DOUBLE)
	call salloc (onorm, ncoeff ** 2, TY_DOUBLE)

	k1 = CV_RANGE(cv)
	k2 = k1 * CV_MAXMIN(cv)

	# Fill normalization matrix after zeroing it. 
	call amovkd (0.0d0, Memd[norm], ncoeff ** 2)
	do n = 1, ncoeff {
	    k = n - 1
	    do i = 0, k {
		bc = dcv_bcoeff (k, i)
		index = norm + i * ncoeff + k
		Memd[index] =  bc * (k1 ** i) * (k2 ** (k-i))
	    }
	}

	# Multiply the input transformation matrix by the normalization
	# matrix.
	call cv_mmuld (Memd[norm], elm, Memd[onorm], ncoeff)
	call amovd (Memd[onorm], elm, ncoeff ** 2)

	call sfree (sp)
end


# CV_BCOEFF -- calculate and return binomial coefficient as function value.

double procedure dcv_bcoeff (n, i)

int	n
int	i

double	dcv_factorial()

begin
	if (i == 0)
	    return (1.0d0)
	else if (n == i)
	    return (1.0d0)
	else
	    return (dcv_factorial (n) / (dcv_factorial (n - i) *
	        dcv_factorial (i)))
end


# CV_FACTORIAL -- calculate factorial of argument and return as function value.

double procedure dcv_factorial (n)

int	n

int	i
double	fact

begin
	if (n == 0)
	    return (1.0d0)
	else {
	    fact = 1.0d0
	    do i = n, 1, -1
	        fact = fact * double (i)
	    return (fact)
	}
end


# CV_MMUL -- Matrix multiply.

procedure cv_mmuld (a, b, c, ndim)

double   a[ndim,ndim]            #I left input matrix
double   b[ndim,ndim]            #I right input matrix
double   c[ndim,ndim]            #O output matrix
int     ndim                    #I dimensionality of system

int     i, j, k
double   v

begin
        do j = 1, ndim
            do i = 1, ndim {
                v = double(0.0)
                do k = 1, ndim
                    #v = v + a[k,j] * b[i,k]
                    v = v + a[k,j] * b[i,k]
                c[i,j] = v
            }
end


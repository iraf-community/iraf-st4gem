include	<math/curfit.h>

include	"dcurfitdef.h"

# CVPSERR -- Convert errors of Legendre or Chebyshev coeffecients to errors
# of power series.

# Jun-27-1994  J.-C. Hsu  adapted from dcvpower

procedure cvpserr (cv, ps_err, ncoeff)

pointer	cv		# Pointer to curfit structure
int	ncoeff		# Number of coefficients in fit
double	ps_err[ARB]	# input/output: orthogonal polynomial coefficients 
			# errors (input); power series coefficients errors 
			# (output)

pointer	sp, elm
int	function
int	dcvstati()

begin
	function = dcvstati (cv, CVTYPE)

	if (function != LEGENDRE && function != CHEBYSHEV) {
	    call eprintf (
		"Cannot convert coefficient errors - wrong function type\n")
	    return
	}

	call smark (sp)
	call salloc (elm, ncoeff ** 2, TY_DOUBLE)

	call amovkd (0.0d0, Memd[elm], ncoeff ** 2)

	switch (function){
	case (LEGENDRE):
	    call legen_err (Memd[elm], ps_err, ncoeff)
	case (CHEBYSHEV):
	    call cheby_err (Memd[elm], ps_err, ncoeff)
	}

	# Normalize 
	call norm_err (cv, ps_err, ncoeff)

	call sfree (sp)
end


# LEGEN_ERR -- Convert errors of Legendre coeffecients to errors of 
# power series coefficients.
# Adapted from cv_legen

procedure legen_err (matrix, ps_err, ncoeff)

int	ncoeff
double	matrix[ncoeff, ncoeff]
double	ps_err[ARB]

int	s, n, r, i
double	sigsq

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

	# Multiply matrix columns by curfit coefficients variance and sum.
	do n = 1, ncoeff {
	    sigsq = 0.0
	    do i = 1, ncoeff
	        sigsq = sigsq + (matrix[i,n] * ps_err[i])**2
	    ps_err[n] = sqrt(sigsq)
	}
end

# CHEBY_ERR -- Convert errors of Chebyshev coeffecients to errors of 
# power series coefficients.
# Adapted from cv_cheby

procedure cheby_err (matrix, ps_err, ncoeff)

int	ncoeff				# Number of coefficients
double	matrix[ncoeff, ncoeff]		# Work array for matrix elements
double	ps_err[ARB]			

int	s, n, m, i
double	sigsq

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

	# Multiply matrix columns by curfit coefficients variance and sum.
	do n = 1, ncoeff {
	    sigsq = 0.0
	    do i = 1, ncoeff
	        sigsq = sigsq + (matrix[i,n] * ps_err[i])**2
	    ps_err[n] = sqrt(sigsq)
	}
end


# NORM_ERR -- Return coefficients errors scaled to full data range.
# Adpated from cv_normalize

procedure norm_err (cv, ps_err, ncoeff)

pointer	cv			# Pointer to curfit structure
int	ncoeff			# Number of coefficients in fit
double	ps_err[ARB]	

pointer	sp, elm, index
int	n, i, k
double	k1, k2, bc, sigsq

double	dcv_bcoeff()

begin
	call smark (sp)
	call salloc (elm, ncoeff ** 2, TY_DOUBLE)

	k1 = CV_RANGE(cv)
	k2 = k1 * CV_MAXMIN(cv)

	# Fill matrix, after zeroing it 
	call amovkd (0.0d0, Memd[elm], ncoeff ** 2)
	do n = 1, ncoeff {
	    k = n - 1
	    do i = 0, k {
		bc = dcv_bcoeff (k, i)
		index = elm + k * ncoeff + i
		Memd[index] = (bc * ps_err[n] * (k1 ** i) * (k2 ** (k-i)))**2
	    }
	}

	# Now sum along matrix columns to get coefficient of individual 
	# powers of x.
	do n = 1, ncoeff {
	   sigsq = 0.0d0
	   do i = 1, ncoeff {
	       index = elm + (n-1) + (i-1) * ncoeff
	       sigsq = sigsq + Memd[index]
	    }
	    ps_err[n] = sqrt(sigsq)
	}

	call sfree (sp)
end

###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	call geterrmat(nfree, fpar, errmat)
#
#  Description:	A procedure to get the diagonal elements of the error matrix
#		for estimating 1 sigma errors on the parameters.
#
#  Arguments:	int	nfree	- number of free parameters
#		real	fpar[ARB]	- array of free parameters
#
#  Returns:	real	errmat[MAXFREE,MAXFREE]	-  the error matrix
#
#  Notes:	Information shared in common blocks defined in "specfit.com".
#
#  History:	June 1989	Gerard Kriss
#		11/1/89	gak	Added automated optimization of step size
#		11/2/89	gak	Revised mixed partials calculation in curvmat
#		11/3/89	gak	Correction for dchi=0 in step size calculation
#
###########################################################################

include "specfit.h"

procedure geterrmat(nfree, fpar, errmat)
int	nfree
real	fpar[ARB]
real	errmat[MAXFREE,MAXFREE]

int	j, k
real	det
double	alpha[MAXFREE, MAXFREE], beta[MAXFREE], delta[MAXFREE]

include "specfit.com"

begin

# Evaluate the curvature matrix.  Curvmat will not divide by the delta's used
#	in calculating the derivatives to avoid FP overflows and underflows.
#	The delta's must be multiplied in explicitly in this routine.
	call curvmat(nfree, fpar, alpha, beta, delta)

# Invert the curvature matrix to get the error matrix.
	call matinv(alpha, nfree, det)

# Fill the elements of the error matrix.
	for ( j = 1; j <= nfree; j = j + 1 ) {
	    for ( k = 1; k <= nfree; k = k + 1 ) {
		errmat[j,k] = (delta[j] * (alpha[j,k]) ) * delta[k]
	    }
	}

end

# CURVMAT - routine to calculate the curvature matrix.
#	Based in part on Bevington's routine CHIFIT (p. 229).
#
procedure curvmat(nfree, fpar, alpha, beta, delta)
int	nfree
real	fpar[ARB]
double	alpha[MAXFREE, MAXFREE]
double	beta[ARB]
double	delta[ARB]

int	j, k
real	aj, ak
real	chisq1, chisq2, chisq3, chisq4, chisq5, dchi

include "specfit.com"

begin

# Fill array of initial step sizes for derivative evaluation.
	for ( j = 1; j <= nfree; j = j + 1) {
		delta[j] = step[ iptr[j] ] / 2.
	}

# Calculate initial chi-square
	call chispec(nfree, fpar, chisq1)

# Determine optimum step sizes to produce delta chisquares of unity.
# 11/1/89 - Approximate as an independent variable in a parabola
	for ( j = 1; j <= nfree; j = j + 1) {
		aj = fpar[j]
		fpar[j] = aj + delta[j]
		call chispec(nfree, fpar, chisq2)
		fpar[j] = aj - delta[j]
		call chispec(nfree, fpar, chisq3)
		dchi = abs(chisq2 + chisq3 - 2. * chisq1)
		if ( dchi > 4.e-4 ) {
			delta[j] = delta[j] * 0.5 / sqrt( dchi )
		} else {
			if ( debug ) {
				call eprintf(" j = %d, dchi = %10.4f\n")
				call pargi(j)
				call pargr(dchi)
			}
			delta[j] = delta[j] * 20.
		}
		fpar[j] = aj
	}

#Generate the curvature matrix (Based on Bevington's CHIFIT, p. 229)
# 11/2/89 Revise calculation of mixed partials to be symmetric about minimum.
#	  This is computationally more expensive, but a better approximation.
	for ( j = 1; j <= nfree; j = j + 1) {
		aj = fpar[j]
		fpar[j] = aj + delta[j]
		call chispec(nfree, fpar, chisq2)
		fpar[j] = aj - delta[j]
		call chispec(nfree, fpar, chisq3)
		fpar[j] = aj
		alpha[j,j] = (chisq2 - 2.*chisq1 + chisq3) / 2.
		beta[j] = - (chisq2 - chisq3) / 4.
		if ( debug ) {
			  call printf("j=%2d %10.4f %10.4f %10.4f delta[j]=%10.3e alpha[j,j]=%10.3e\n")
				call pargi(j)
				call pargr(chisq1)
				call pargr(chisq2)
				call pargr(chisq3)
				call pargd(delta[j])
				call pargd(alpha[j,j])
		}
		if ( alpha[j,j] < 0. ) {
			if ( -alpha[j,j]/chisq1 < 4.e-7 ) {
				alpha[j,j] = -alpha[j,j]
			}
		}
		for ( k = j + 1; k <= nfree; k = k + 1) {
			aj = fpar[j]
			ak = fpar[k]
			fpar[j] = aj + delta[j] / 2.
			fpar[k] = ak + delta[k] / 2.
			call chispec(nfree, fpar, chisq2)
			fpar[k] = fpar[k] - delta[k]
			call chispec(nfree, fpar, chisq3)
			fpar[j] = fpar[j] - delta[j]
			call chispec(nfree, fpar, chisq4)
			fpar[k] = fpar[k] + delta[k]
			call chispec(nfree, fpar, chisq5)
			fpar[k] = ak
			fpar[j] = aj

			alpha[k,j] = (chisq2 - chisq3 - chisq5 + chisq4) / 2.
	#		if ( abs(alpha[k,j])/chisq1 < 4.e-7 ) {
	#			if ( debug ) {
	#				call printf("k,j=%3d, %3d.  Small alpha[k,j] = %14.8e\n")
	#				call pargi(k)
	#				call pargi(j)
	#				call pargd(alpha[k,j])
	#			}
	#			alpha[k,j] = 0.
	#		}
			if ( k == 1 && debug ) {
			  call printf("k,j=%2d,%2d alpha[k,j] = %14.8e\n")
				call pargi(k)
				call pargi(j)
				call pargd(alpha[k,j])
			}
			alpha[j,k] = alpha[k,j]
			if ( debug ) {
			   call printf("k=%3d %10.4f %10.4f %10.4f %10.4f delta[k]=%10.3e alpha[j,k]=%10.3e\n")
					call pargi(k)
					call pargr(chisq2)
					call pargr(chisq3)
					call pargr(chisq4)
					call pargr(chisq5)
					call pargd(delta[k])
					call pargd(alpha[j,k])
			}
		}
	}

# Make sure parameters are reset to their original values in the common block
	call update(nfree, fpar)
end

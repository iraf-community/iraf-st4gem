###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	specfit
#
#  Description:	SPECFIT is an IRAF task for fitting complex continua with
#		mulitcomponents including emission lines.
#
#  Arguments:	int	nfree
#		real	fpar[ARB], fstep[ARB], maxerr[ARB]
#		real	tol
#		int	maxiter
#		extern	fx
#		real	fx
#		real	chisq
#
#  Returns:	none
#
#  Notes:	fx() is an external user-supplied function
#
#
#  History:	May 1989	Gerard Kriss
#				Wrote major portions in SPP, based on the routine
#				CURFIT (Bevington 1969, Data Reduction and Error
#				Analysis for the Physical Sciences, pp. 237-239.)
#
#
###########################################################################

include "specfit.h"

real procedure marquadt(nfree, fpar, fstep, maxerr, tol, maxiter, fx, chisq)
int	nfree
real	fpar[ARB], fstep[ARB], maxerr[ARB]
real	tol
int	maxiter
extern	fx
real	chisq

int	j, k, niter, miter, curve, done
real	chisq1, chisq2
real	det, error, newpar[MAXFREE]
double	flambda
double	alpha[MAXFREE, MAXFREE], aprime[MAXFREE, MAXFREE]
double	beta[MAXFREE], delta[MAXFREE]

bool	err_from_model
bool	interact
bool	debug
int	nlogfd
int	logfd[5]
common /users/ err_from_model, interact, debug, nlogfd, logfd



begin

	flambda = 1.0
	call fx(nfree, fpar, chisq1)        # Get initial chi-square
	niter = 0
	done = 0
	curve = 1
	while ( (niter < maxiter) && (done == 0) ) {
	niter = niter + 1

# Before starting, impose user-defined limits on parameters.
	call setlim(nfree, fpar)

# Evaluate curvature matrix, alpha, and gradient vector, beta.
# Remember, a delta[j] * delta[k] is missing from the denominator of
# each alpha[j,k].  Multiply it in later.
	call curvmat(nfree, fpar, alpha, beta, delta)

	chisq2 = chisq1 + 1.e-3
	miter = 0
	while ( (chisq2 > chisq1) && (miter < 2) ) {
	miter = miter + 1
	flambda = flambda * 10.
# Add Marquadt parameter to the diagonal of alpha and check that curvature
# is positive.
	curve = 1
	for ( j = 1; j <= nfree; j = j+ 1) {
		for ( k = 1; k <= nfree; k = k + 1) {
			aprime[j,k] = alpha[j,k]
		}
		if ( alpha[j,j] < 0. ) {
			#alpha[j,j] = - alpha[j,j]
			curve = -1
		}
		aprime[j,j] = alpha[j,j] * (1. + flambda)
	}
	call matinv(aprime, nfree, det)

# Compute changes to current free parameters
	for ( j = 1; j <= nfree; j = j + 1) {
	    newpar[j] = 0.
	    for ( k = 1; k <= nfree; k = k + 1) {
		#call printf("k=%2d beta[k]=%g aprime[j,k]=%g delta[k]=%g.\n")
		#call pargi(k)
		#call pargr(beta[k])
		#call pargr(aprime[j,k])
		#call pargr(delta[k])
	        newpar[j] = newpar[j] +  (beta[k] * aprime[j,k])
	    }
		#call printf("j=%2d newpar[j]=%g delta[j]=%g.\n")
		#call pargi(j)
		#call pargr(newpar[j])
		#call pargr(delta[j])
	    newpar[j] = newpar[j] * delta[j] + fpar[j]
	}
	call fx(nfree, newpar, chisq2)
	if ( debug ) {
		call printf("niter = %3d.  chisq1, chisq2 = %g %g.\n")
		call pargi(niter)
		call pargr(chisq1)
		call pargr(chisq2)
	}
	}

	if ( chisq2 < chisq1 ) {
	# Chi-square decreased, so decrease flambda and check for done.
		flambda = flambda / 100.
		for ( j = 1; j <= nfree; j = j + 1) {
			fpar[j] = newpar[j]
		}
		done = 1
		if ( curve < 0 ) {
			done = 0
			if ( debug ) {
				call printf("Curvature < 0.\n")
			}
		}
		if ( abs((chisq1 - chisq2)/chisq1) > tol )
			done = 0
		for ( j = 1; j <= nfree && done == 1; j = j + 1) {
			error = (newpar[j] - fpar[j]) / fpar[j]
			if ( abs(error) > maxerr[j] )
				done = 0
		}
		chisq1 = chisq2
	}
	call printf("%d %15.7g\n")	# Report on progress
		call pargi(niter)
		call pargr(chisq1)

	}
	chisq = chisq1
end

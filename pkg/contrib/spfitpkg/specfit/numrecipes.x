###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	specfit
#
#  Description:	mr_solve is an optimization algorithmn that finds the best fit
#		between data and a function ( in this case fspec ) 
#
#  Arguments:	int	nfree
#		real	fpar[ARB], fstep[ARB], maxerr[ARB]
#		real	tol
#		int	maxiter
#		real	chisq
#
#  Returns:	none
#
#
#  History:	Dec 1994	J Grimes	Based on version from Numerical
#						Recipes that is part of iraf
#						libraries
#
#
###########################################################################

include <math.h>
include <mach.h>

include "specfit.h"


real procedure mr_solve(nfree, fpar, fstep, maxerr, tol, maxiter, chisq)
int	nfree
real	fpar[ARB], fstep[ARB], maxerr[ARB]
real	tol
int	maxiter
real	chisq


pointer	Beta,ftry,Da,Covar,Alpha

int 	j,count
real 	mr
int 	iter
int 	done
real	OChisq
real	error

include "specfit.com"

int mod()

begin

iter = 0
done = NO
count = -1

call malloc(Beta,nfree,TY_REAL)
call malloc(ftry,nfree,TY_REAL)
call malloc(Da,nfree,TY_REAL)
call malloc(Covar,nfree*nfree,TY_REAL)
call malloc(Alpha,nfree*nfree,TY_REAL)

# Initialization
mr = 0.001

do j=1, nfree*nfree {
	Memr[Alpha+j-1] = 0.0
}

call mr_eval( fpar,nfree,Memr[Alpha],Memr[Beta],Chisq, fstep )

OChisq = Chisq
call amovr(fpar,Memr[ftry],nfree)

# Main minimization loop
while ( done != YES && iter < maxiter ) {

	iter = iter + 1
	
	call setlim(nfree,fpar)

	if ( mod(iter,4) == 0 ) {
		mr = 0.001
 	}

	call amovr(Memr[Alpha],Memr[Covar],nfree * nfree)
	call amovr(Memr[Beta],Memr[Da],nfree)
	do j = 1, nfree {
		Memr[Covar+(j-1)*(nfree+1)] = Memr[Alpha+(j-1)*(nfree+1)] *
						(1.0 + mr)

	}

	call mat_inverse(nfree,Memr[Covar], Memr[Da]) 

	do j = 1, nfree {
		Memr[ftry+j-1] = Memr[ftry+j-1] + Memr[Da+j-1]
	}

	call mr_eval( Memr[ftry],nfree, Memr[Covar], Memr[Da], Chisq, fstep )

	if ( Chisq < OChisq ) {
		mr = 0.1 * mr

		done = YES
		if ( (OChisq - Chisq)/OChisq > tol ) {
			done = NO
		}

		for ( j = 1; j <= nfree && done == 1; j = j + 1) {
			error = (Memr[ftry+j-1] - fpar[j]) / fpar[j]
			if ( abs(error) > maxerr[j] )
				done = NO
		}

		OChisq = Chisq

		call amovr(Memr[Covar],Memr[Alpha],nfree*nfree)
		call amovr(Memr[Da],Memr[Beta],nfree)
		call amovr(Memr[ftry],fpar,nfree)

	} else {
		mr = 10.0 * mr
		Chisq = OChisq
	}

	call printf("\t\t\t\t%d %15.7g\n")	# Report on progress
		call pargi(iter)
		call pargr(chisq)
	
  } #End MAin While Loop

call mfree(Covar,TY_REAL)
call mfree(Beta,TY_REAL)
call mfree(Alpha,TY_REAL)
call mfree(ftry,TY_REAL)
call mfree(Da,TY_REAL)

end


procedure mr_eval( fpar, nfree, Alpha, Beta, Chisq , fstep )
real 	fpar[MAXFREE]
int 	nfree
real	Alpha[nfree,nfree], Beta[nfree]
real 	Chisq
real 	fstep[ARB]

int 	i,j,k
real 	Wt, Sig2I, Dy
real 	DyDa[MAXFREE]
real 	val1, val2, savej, term


include "specfit.com"

begin	

	do j = 1, nfree {
		do k = 1, j {
			Alpha[j,k] = 0.0
		}
		Beta[j] = 0.0
	}

	chisq = 0.0

	do i = 1, npts {

		if ( Memi[infit+i-1] == 1 ) {

			# Calculate partial derivatives wrt the parameters
			do j = 1, nfree {
				savej = fpar[j]
				fpar[j] = savej + fstep[j]
				call update(nfree,fpar)
				call fspec( Memr[lambdas+i-1], val1 )
				fpar[j] = savej - fstep[j]
				call update(nfree,fpar)
				call fspec(Memr[lambdas+i-1], val2)
				dyda[j] = ( val1 - val2 ) / ( 2. * fstep[j] )
				fpar[j] = savej
			}
			
			# Get fitted weights and calculate Chi-square

			call update(nfree,fpar)
			call fspec(Memr[lambdas+i-1], Memr[fitsp+i-1])	

			if ( err_from_model ) {
			    # Assume Poisson statistics using expected
			    # variance from model. Uses PROS approximation
			    # sigma' = 1 + sqrt(n+0.75), but with n
			    # from the model, not the data
				term = (Memr[spectrum+i-1] - Memr[fitsp+i-1])
				if ( Memr[fitsp+i-1] <= -0.75 ) {
				    Sig2I = 3.48205
				} else {
				    Sig2I = (1.+sqrt(Memr[fitsp+i-1]+0.75))**2
		    		    term = term*term / Sig2I
				}
				chisq = chisq + term
			} else {
			    # Assume Gaussian statistics with errors
			    # supplied by the user
				term = (Memr[spectrum+i-1] - Memr[fitsp+i-1]) / Memr[errors+i-1]
				Sig2I = Memr[errors+i-1] * Memr[errors+i-1]
				chisq = chisq + term**2
			}
		
			# Calculate elements of alpha and beta matrices
			Dy = Memr[spectrum+i-1] - Memr[fitsp+i-1]
			do j = 1, nfree {
				Wt = DyDa[j] / Sig2I
				do k = 1, j {
 				  Alpha[j,k] = Alpha[j,k] + Wt * DyDa[k]
				}
				Beta[j] = Beta[j] + Dy * Wt 
			}

		} else {
			Memr[fitsp+i-1] = 0.
		}
	}

	do j = 2, nfree {
		do k = 1, j-1 {
			Alpha[k,j] = Alpha[j,k]
		}
	}

end

include "pltrans.h"

# SERIESFIT -- Fits a series of polynomials to a set data: the specified
#		degree and lesser degrees.  Returns chisqr for each degree
#		for the user to determine best fit.

procedure seriesfit(xdat, ydat, ysig, npt, nterm, mode, coeff, chisqr)

real	xdat[ARB]		# i: xdata array
real	ydat[ARB]		# i: ydata array
real	ysig[ARB]		# i: y data errors
int	npt			# i: number of data points
int	nterm			# i: number of terms in polynomial (degree+1)
int	mode			# i: fitting mode
real	coeff[MAXFIT,ARB]	# o: array of coeffecients for fits
real	chisqr[ARB]		# o: Chisqr values for fits

int	ic, jc, minterm
pointer	cf
real	chi

begin

	# Perform up to MAXFIT polynomial fits of nterm terms or less

	minterm = max( 1, nterm-MAXFIT )
	call malloc( cf, nterm, TY_REAL)

	do ic = minterm, nterm {
	   
	   call polfit( xdat, ydat, ysig, npt, ic, mode, Memr[cf], chi )

	   do jc = 1, nterm
	      coeff[ic,jc] = Memr[cf+jc-1]

	   chisqr[ic] = chi

	}

	call mfree( cf, TY_REAL)
end

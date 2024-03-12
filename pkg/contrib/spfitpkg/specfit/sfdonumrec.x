###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	procedure sfdonumrec(nfree, fpar, chisq)
#
#  Description:	SFDONUMREC performs a fit using the current best guesses for
#		the parameters with a simplex routine and returns the new chisq
#
#  Arguments:	int	nfree		- Number of free params
#		real	fpar[ARB]	- Array containing the free params
#
#  Returns:	real	chisq		- Chisquare for the current fit
#
#  Notes:	Shares data in "specfit.com"
#
#  History:	July 	1994	Grimes  Taken from Numerical Recipes in
#					Iraf libraries
#
###########################################################################

include	"specfit.h"

procedure sfdonumrec(nfree, fpar, chisq)
int	nfree
real	fpar[ARB]
real	chisq

int	i
real	fstep[MAXFREE], ftol[MAXFREE]


include	"specfit.com"

begin
#call printf("Entered sfdonumrec.\n")
	call freezepar(nfree, fpar)
	if ( nfree > 0 ) {
		for ( i = 1; i <= nfree; i = i + 1 ) {
			fstep[i] = step[ iptr[i] ]
			ftol[i] = ptol[ iptr[i] ]
#call printf("%2d %2d %g %g\n")
 # call pargi(i)
 # call pargi(iptr[i])
 # call pargr(fstep[i])
 # call pargr(ftol[i])
		}
		call mr_solve(nfree,fpar,fstep,ftol,tolerance,itr,chisq)
		call update(nfree, fpar)
	} else {
		call chispec(nfree, fpar, chisq)
	}
end






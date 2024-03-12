###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University

#  Synopsis:	procedure sfdosimplex(nfree, fpar, chisq)

#  Description:	SFDOSIMPLEX performs a fit using the current best guesses for
#		the parameters with a simplex algorithm
#		and returns chisquare and the updated parameters

#  Arguments:	int	nfree		- Number of free params
#		real	fpar[ARB]	- Array containing the free params

#  Returns:	real	chisq		- Chisquare for the current fit

#  Notes:	Shares data in "specfit.com"

#  History:	May	1989	Gerard Kriss

###########################################################################

include	"specfit.h"

procedure sfdosimplex(nfree, fpar, chisq)
int	nfree
real	fpar[ARB]
real	chisq

int	i
real	fstep[MAXFREE], ftol[MAXFREE]

extern chispec

include	"specfit.com"

begin
#call printf("Entered sfdosimplex.\n")
	call freezepar(nfree, fpar)
	if ( nfree > 0 ) {
		for ( i = 1; i <= nfree; i = i + 1 ) {
			fstep[i] = step[ iptr[i] ]
			ftol[i] = ptol[ iptr[i] ]
#call printf("%2d %2d %g %g\n")
#  call pargi(i)
#  call pargi(iptr[i])
#  call pargr(fstep[i])
#  call pargr(ftol[i])
		}
		call simplex(nfree,fpar,fstep,ftol,tolerance,itr,chispec,chisq)
		call update(nfree, fpar)
	} else {
		call chispec(nfree, fpar, chisq)
	}
end

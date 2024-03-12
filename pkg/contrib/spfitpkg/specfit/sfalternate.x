###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	procedure sfalternate(nfree, fpar, chisq)
#
#  Description:	SFALTERNATE performs a fit using the current best guesses for
#		the parameters alternating between 3 iterations with a Simplex
#		algorithm and 3 with a Marquadt minimization algorithm.
#		Returns chisquare and the updated parameters
#
#  Arguments:	int	nfree		- Number of free params
#		real	fpar[ARB]	- Array containing the free params
#
#  Returns:	real	chisq		- Chisquare for the current fit
#
#  Notes:	Shares data in "specfit.com"
#
#  History:	May	1989	Gerard Kriss

###########################################################################

include	"specfit.h"

procedure sfalternate(nfree, fpar, chisq)
int	nfree
real	fpar[ARB]
real	chisq

int	niter, maxiter_save

extern chispec

include	"specfit.com"

begin
	maxiter_save = itr
	niter = 0
	if ( maxiter_save == 0 ) {
		call chispec(nfree, fpar, chisq)
	}
	while ( niter < maxiter_save ) {
		niter = niter + 1
		itr = 3
		call sfdosimplex(nfree, fpar, chisq)
		itr = 5
		call sfdonumrec(nfree, fpar, chisq)
	}
	itr = maxiter_save
end

###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	call chispec(nfree, fpar, chisq)
#
#  Description:	A procedure to evaluate chi-square for fits of
#		multi-component spectra for SPECFIT
#
#  Arguments:	int	nfree	- number of free parameters
#		real	fpar	- array of free parameters
#
#  Returns:	real	chisq	- returned value of chisquare
#
#  Notes:	Information shared in common blocks defined in "specfit.com".
#		Spectral data obtained from /spec/
#
#  History:	May 1989	Gerard Kriss
#               Feb 1992	H. Ferguson - added option to use expected
#				variance from model in chisq calculation
#               May 1995        H. Ferguson - Modified to use 
#                               sigma=1+sqrt(n+0.75) when variance from model
#
###########################################################################

include "specfit.h"

procedure chispec(nfree, fpar, chisq)
int	nfree
real	fpar[ARB]
real	chisq

int	i
real	term

include "specfit.com"

begin
#call printf("Entered chispec.\n")

# Update free parameters in the common block
	call update(nfree, fpar)
#call printf("Updated params.\n")

	chisq = 0.
	# if ( debug ) {
	# 	call printf("chisq: errors from model: %g %g\n")
	# 	call pargb(err_from_model)
	# 	call pargb(debug)
	# }
	if ( err_from_model ) {   
	    # Assume Poisson statistics using expected variance from model
	    # Uses PROS approximation sigma' = 1 + sqrt(n+0.75), but with n
            # given by the model, not the data.
	    for ( i = 0; i < npts; i = i + 1) {
		if ( Memi[infit+i] == 1 ) {
		    call fspec(Memr[lambdas+i], Memr[fitsp+i])
		    term = (Memr[spectrum+i] - Memr[fitsp+i])
		    term = term*term / (1. + sqrt(Memr[fitsp+i]+0.75))**2
		    chisq = chisq + term
		} else {
			Memr[fitsp+i] = 0.
		}
	    }
	} else {
	    # Assume Gaussian statistics with errors supplied by the user
	    for ( i = 0; i < npts; i = i + 1) {
		if ( Memi[infit+i] == 1 ) {
		    call fspec(Memr[lambdas+i], Memr[fitsp+i])
		    term = (Memr[spectrum+i] - Memr[fitsp+i]) / Memr[errors+i]
		    chisq = chisq + term**2
		} else {
			Memr[fitsp+i] = 0.
		}
	    }
        }
	if ( debug ) {
		call printf("chisq = %g\n")
			call pargr(chisq)
	}
end

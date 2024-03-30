# ASPECFUNK -- Evaluate chisqr for the given parameter set.  Used in 
# aspecfit

real procedure specfunk( par )

real	par[ARB]
#--
int	npts
real	chi2
char	model[SZ_LINE]

# Include common block for communication with ASPECFIT
include	"../lib/amoebafit.h"

begin
	# Copy model string so it can be updated
	call strcpy( fitmodel, model, SZ_LINE)
	call insertpar( par, model )

	# Evaluate chisqr
	call specchi2( model, fitdata, chi2, npts )

	return( chi2 )
end

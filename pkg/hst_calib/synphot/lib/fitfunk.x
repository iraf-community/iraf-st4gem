# FITFUNK -- Evaluate chisqr for the given parameter set.  Used in 
# Amoebafit

real procedure fitfunk( par )

real	par[ARB]
#--
int	npts
real	chi2
char	model[SZ_LINE]

# Include common block for communication with AMOEBAFIT
include	"amoebafit.h"

begin
	# Copy model string so it can be updated
	call strcpy( fitmodel, model, SZ_LINE)
	call insertpar( par, model )

	# Evaluate chisqr
	call evalchi2( model, fitdata, chi2, npts )

	return( chi2 )
end

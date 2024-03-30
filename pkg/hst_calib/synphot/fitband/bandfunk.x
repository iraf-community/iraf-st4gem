include "../plspec/plspec.h"

# BANDFUNK -- Evaluate chisqr for the given parameter set.  Used in 
# Amoebafit.  Fitband version

real procedure bandfunk( par )

real	par[ARB]
#--
int	npts
real	chi2, band[MAXPHOT]
char	model[SZ_LINE]

# Include common block for communication with AMOEBAFIT
include	"../lib/amoebafit.h"

begin
	# Copy model string so it can be updated
	call strcpy( fitmodel, model, SZ_LINE)
	call insertpar( par, model )

	# Evaluate chisqr
	call bandchi2( model, fitdata, chi2, band, npts )

	return( chi2 )
end

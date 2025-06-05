include <iraf77.h>

# UGPGP[IR]-- Set a GIO parameter.
# GIO parameters are of type integer or real.
# The only character string GIO parameter is xtickformat, and we do not 
# support it!

procedure ugpgpr (gp, param, value, istat)

pointer	gp			# graphics descriptor
int	param			# parameter to be gotten
real	value			# value of parameter
int	istat

begin
	istat = ER_OK

	iferr (call gsetr (gp, param, value))
	    istat = ER_GRAPHSETPAR
end

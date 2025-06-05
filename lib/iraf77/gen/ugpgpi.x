include <iraf77.h>

# UGPGP[IR]-- Set a GIO parameter.
# GIO parameters are of type integer or real.
# The only character string GIO parameter is xtickformat, and we do not 
# support it!

procedure ugpgpi (gp, param, value, istat)

pointer	gp			# graphics descriptor
int	param			# parameter to be gotten
int	value			# value of parameter
int	istat

begin
	istat = ER_OK

	iferr (call gseti (gp, param, value))
	    istat = ER_GRAPHSETPAR
end

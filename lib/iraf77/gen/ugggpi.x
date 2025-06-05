include <iraf77.h>

# UGGGP[IR]-- Get a GIO parameter.
# GIO parameters are of type integer or real.
# The only character string GIO parameter is xtickformat, and we do not 
# support it!

procedure ugggpi (gp, param, value, istat)

pointer	gp			# graphics descriptor
int	param			# parameter to be gotten
int	value			# real value of parameter
int	istat

int	gstati()		# function to return GIO parameter

begin
	istat = ER_OK

	iferr (value = gstati (gp, param))
	    istat = ER_GRAPHGETPAR
end

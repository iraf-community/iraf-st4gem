include <iraf77.h>

# UGRSET Reset internal GIO parameters to their default values by class 
# of parameters rather than one by one with ugset

procedure ugrset (gp, parcls, istat)

pointer	gp			# graphics descriptor
int	parcls			# parameter class to reset
int	istat

begin
	iferr {
	    call greset (gp, parcls)
	} then
	    istat = ER_GRAPHILLCODE	# reset code is illegal
	else
	    istat = ER_OK
end

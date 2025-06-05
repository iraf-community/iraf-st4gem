include <iraf77.h>

# UGCLOS  -- Close a graphics station

procedure ugclos (gp, istat)

pointer	gp
int	istat

begin
	istat = ER_OK

	# Close graphics station.
	iferr (call gclose (gp))
	   istat = ER_GRAPHCLOSE
end

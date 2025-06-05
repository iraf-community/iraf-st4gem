include <iraf77.h>

# UGCLRS  -- Clear a graphics station
#
# If the device is a CRT, the screen is cleared.  If hardcopy, ugclrs() 
# advances the frame (starts a new page).  

procedure ugclrs (grdscr, istat)

pointer	grdscr		# Graphics descriptor
int	istat

begin
	istat = ER_OK

	# Clear graphics station.
	iferr (call gframe (grdscr))
	   istat = ER_GRAPHCLEAR
end

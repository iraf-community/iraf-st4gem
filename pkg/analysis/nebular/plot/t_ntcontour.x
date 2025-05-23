# Copyright restrictions apply - see st4gem$copyright.stsdas

include <gset.h>
include	"../neberr.h"
include	"ntplot.h" 

define	DEBUG	false

#---------------------------------------------------------------------6 Jul 97--
.help t_ntcontour.x Mar97 nebular
.ih
NAME
t_ntcontour - Task to generate contours of temperature/density sensitive ratios
.endhelp
#-------------------------------------------------------------------------------
# T_NTCONTOUR - Task to generate contours of temperature/density sensitive 
#		ratios.  

procedure t_ntcontour ()

#  Declarations
bool	clgetb()		# fetch CL parameter TY_BOOL
char	device[SZ_FNAME]	# graphics device
pointer gopen()                 # Returns pointer to graphics structure
pointer	gp			# graphics descriptor
bool	interactive		# interactive session?
pointer	pl			# plot state descriptor
pointer	plt_alloc()		# allocate plot descriptor

errchk	get_cl_pltpars, gopen, nt_contour, plot_setup

begin
#	if (DEBUG)
#	    call memlog ("---------- Start ntcontour ----------")

	# Fetch initial plot attributes from the CL.
	pl = plt_alloc (INDEFI)
	call get_cl_pltpars (pl)

	# Set up graphics descriptor.
	call clgstr ("device", device, SZ_FNAME)
	PL_APPEND(pl) = clgetb ("append")
	if (PL_APPEND(pl))
	    gp = gopen (device, APPEND, STDGRAPH)
	else
	    gp = gopen (device, NEW_FILE, STDGRAPH)

	PL_GP(pl) = gp
	call gseti (gp, G_NTITLELINES, 3)

	# Plot curves. 
	interactive = clgetb ("interactive")
	call nt_contour (pl, interactive)
	call plt_free (pl)

#	if (DEBUG)
#	    call memlog ("---------- End ntcontour ----------")
end



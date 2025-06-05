include <iraf77.h>

# UGOPEN  -- Open a graphics workstation
#
#  Ugopen() is restricted to open the STDGRAPH stream only.  The device
#  name is passed as a Fortran 77 character string and must resolve to an
#  entry in the graphcap file.  If the device is specified as "stdgraph", 
#  "stdplot", "stdimage", or "stdvdm", the value of the environment parameter 
#  of the same name will be used.  Otherwise, the string is assumed to be a 
#  recognized device.

procedure ugopen (mode, f77dev, gp, istat)

int	mode		# Graphics stream open mode (APPEND or NEWFILE)
%      character*(*)	f77dev
pointer	gp		# Pointer to graphics
int	istat		# Status

int	txtlen
pointer	sp, device
int	fd		# Graphics stream descriptor


pointer	gopen()

begin
	istat = ER_OK

	# Check that mode is OK.
	if (mode != NEW_FILE && mode != APPEND) {
	   istat = ER_BADGRAPHMODE
	   return
	}

	# Convert F77 device string to SPP device string
	txtlen = len (f77dev)
	call smark (sp)
	call salloc (device, txtlen, TY_CHAR)
	call f77upk (f77dev, Memc[device], txtlen)

	# Always use STDGRAPH stream to avoid confusion;
	# The user can still use "stdplot", etc. to write to the 
	# appropriate device
	fd = STDGRAPH

	# Open the graphics stream
	iferr (gp = gopen (Memc[device], mode, fd)) {
	   istat = ER_GRAPHOPEN
	}

	call sfree (sp)
end

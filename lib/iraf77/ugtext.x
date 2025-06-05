include <iraf77.h>
include	<gset.h>

# UGTEXT -- Draw text.  All textual output is via this primitive.
# Textual output is not subjected to clipping Clipping may be performed at 
# the kernel level if a workstation viewport is defined.  

procedure ugtext (gp, x, y, f77txt, istat)

pointer	gp			# graphics descriptor
real	x, y			# position at which text is to be drawn
%	character*(*) f77txt
int	istat			# return status

pointer	sp, text
int	txtlen			# size of input string

begin
	istat = ER_OK

	# convert text to SPP string
	txtlen = len (f77txt)
	call smark (sp)
	call salloc (text, txtlen, TY_CHAR)
	call f77upk (f77txt, Memc[text], txtlen)

	# Plot text 
	# Since we set all GIO parameters internally we don't have
	# to worry about the format.
	iferr (call gtext (gp, x, y, Memc[text], ""))
	   istat = ER_GRAPHTEXT

	call sfree (sp)
end

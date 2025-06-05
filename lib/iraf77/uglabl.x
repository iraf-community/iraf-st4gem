include <iraf77.h>
include	<gset.h>

# UGLABL -- Draw and label the axes of the plot (normally the viewport
# boundary).  This is done in two steps.  First we compute all the required
# parameters, and then we draw and label the axes.  Up to four axes can be
# drawn.  To simplify matters, all four axes are treated equally and
# independently.  The axes are drawn a tick at a time in world coordinates.

procedure uglabl (grdscr, f77ttl, f77xlb, f77ylb, istat)

pointer	grdscr			# graphics descriptor
%	character*(*) f77ttl
%	character*(*) f77xlb
%	character*(*) f77ylb
int	istat			# return status

int	txtlen
pointer	sp, title, xlabel, ylabel

begin
	istat = ER_OK

	# Convert Fortran 77 strings to SPP strings
	call smark (sp)
	txtlen = len (f77ttl)
	call salloc (title, txtlen, TY_CHAR)
	call f77upk (f77ttl, Memc[title], txtlen)
	if (Memc[title] != EOS)
	    call gseti (grdscr, G_DRAWTITLE, YES)
	else
	    call gseti (grdscr, G_DRAWTITLE, NO)

	txtlen = len (f77xlb)
	call salloc (xlabel, txtlen, TY_CHAR)
	call f77upk (f77xlb, Memc[xlabel], txtlen)
	if (Memc[xlabel] != EOS)
	    call gseti (grdscr, G_XLABELAXIS, YES)
	else
	    call gseti (grdscr, G_XLABELAXIS, NO)

	txtlen = len (f77ylb)
	call salloc (ylabel, txtlen, TY_CHAR)
	call f77upk (f77ylb, Memc[ylabel], txtlen)
	if (Memc[ylabel] != EOS)
	    call gseti (grdscr, G_YLABELAXIS, YES)
	else
	    call gseti (grdscr, G_YLABELAXIS, NO)

	# Do the actual work now
	iferr {
	    call glabax (grdscr, 
		Memc[title], Memc[xlabel], Memc[ylabel])
	} then
	    istat = ER_GRAPHLABEL

	call sfree (sp)
end

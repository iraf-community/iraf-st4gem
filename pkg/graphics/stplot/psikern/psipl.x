include <gki.h>
include <gset.h>
include "psi.h"

define  MAX_LTYPES      3       # max software line type patterns (excl. solid)
define  MAX_LSEGMENTS   4       # max line segments per pattern
define  LT_OFFSET       1       # Offset to be subtracted from ltype code

#---------------------------------------------------------------------------
.help psi_polyline 1May92 plot
.ih
NAME
psi_polyline -- Draw a polyline.  
.ih
USAGE
call psi_polyline (p, npts)
.ih
ARGUMENTS
.ls p (short[ARB])
The points defining the line.
.le
.ls npts (int)
The number of points, i.e., (x,y) pairs in p.
.le
.ih
DESCRIPTION
The polyline is defined by the array of points P, consisting of 
successive (x,y) coordinate pairs.  The first point is not plotted but
rather defines the start of the polyline.  The remaining points define
line segments to be drawn.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_polyline (p, npts)

short   p[ARB]                  # points defining line
int     npts                    # number of points, i.e., (x,y) pairs

# Declarations.
int i       # Generic.
int len_p   # The total number of values in the P array.
int x, y    # The x, y point.

pointer pl  # Pointer to the polyline attribute descriptor.

include "psi.com"

begin
        
        if (g_debug){
            call eprintf ("psi_polyline: Drawing one of %d points.\n")
            call pargi (npts)
        }
        
        # If not enough points, just forget it.
        if (npts < 2)
            return
        
        len_p = npts * 2
        
        # Keep track of the number of drawing instructions since the last frame
        # clear.
        g_ndraw = g_ndraw + 1

	# Set DCS page indicator.
	if (g_bpage) {
	    call psk_page
	    g_bpage = false
	}
	
        # Update polyline attributes if necessary.
        pl = PSI_PLAP(g_kt)
        call psi_ltype (PL_LTYPE(pl))
        call psi_linewidth (PL_WIDTH(pl))
        
        # Set the appropriate color.  In PostScript, erasure is handled by
        # writing in the same color as the background.
        if (PL_LTYPE(pl) == GL_CLEAR)
            call psi_erase
        else
            call psi_color (PL_COLOR(pl))
        
        # Move to the first point.
        x = p[1]
        y = p[2]
        call psi_move (x, y)
        
        # Draw the polyline.
        for (i=3;  i <= len_p;  i=i+2) {
            x = p[i]
            y = p[i+1]
            call psi_draw (x, y)
        }
        
        # Set the drawing commands to the page.
        call psi_stroke
        
end
#---------------------------------------------------------------------------
# End of psi_polyline
#---------------------------------------------------------------------------

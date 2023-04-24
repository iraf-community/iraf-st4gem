include <gki.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_polymark 1May92 plot
.ih
NAME
psi_polymarker -- Draw a polymarker.  
.ih
USAGE
call psi_polymarker (p, npts)
.ih
ARGUMENTS
.ls p (short[ARB])
The points defining the line.
.le
.ls npts (int)
The number of (x,y) pairs in the p vector.
.le
.ih
DESCRIPTION
The polymarker is defined by the array
of points P, consisting of successive (x,y) coordinate pairs.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_polymarker (p, npts)

short   p[ARB]                  # points defining line
int     npts                    # number of points, i.e., (x,y) pairs

# Declarations.
int i, len_p
int x, y, oldx, oldy

pointer pm

include "psi.com"

begin
        
        if (g_debug) {
            call eprintf ("psi_polymarker: Drawing on of %d points.\n")
            call pargi (npts)
        }
        
        if (npts <= 0)
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
        
        # Update polymarker attributes if necessary.
        pm = PSI_PMAP(g_kt)
        call psi_ltype (PM_LTYPE(pm))
        call psi_linewidth (PM_WIDTH(pm))
        call psi_color (PM_COLOR(pm))
        
        # Draw the polymarker.
        oldx = 0; oldy = 0
        for (i=1;  i <= len_p;  i=i+2) {
            x = p[i]; y = p[i+1]
            if (x != oldx || y != oldy) {
                call psi_move (x, y)
                call psi_draw (x, y)
            }
            oldx = x; oldy = y
        }

        # Set the drawing commands to the page.
        call psi_stroke
        
end
#---------------------------------------------------------------------------
# End of psi_polymarker
#---------------------------------------------------------------------------

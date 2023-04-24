include <gset.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_fillarea 1May92 plot
.ih
NAME
psi_fillarea -- Fill a closed area.
.ih
USAGE
call psi_fillarea (p, npts)
.ih
ARGUMENTS
.ls p (short[ARB])
Points defining a line.
.le
.ls npts (int)
The number of (x,y) pairs in the vector p.
.le
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_fillarea (p, npts)

short   p[ARB]                  # points defining line
int     npts                    # number of points, i.e., (x,y) pairs

# Declarations.
int     i             # Generic.
int     len_p         # The total number of values in the P array.
int     old_pl_color  # The original color of polylines.
int     old_pl_ltype  # The original line type of polylines.
int     x, y          # The coordinates of the current point.

pointer fa            # Pointer to the Fill Area attribute descriptor.
pointer pl            # Pointer to the PolyLine attribute descriptor.

include "psi.com"

begin
        if (g_debug) {
            call eprintf ("psi_fillarea: Drawing area of %d points.\n")
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
        fa = PSI_FAAP(g_kt)
        
        # Set the fill style.
        PSI_FILLSTYLE(g_kt) = FA_STYLE(fa)  
        
        # See if color has changed.  One change is that if the fill type is
        # clear, this corresponds to the default background of the screen,
        # usually white.
        if (FA_STYLE(fa) == GF_CLEAR)
            call psi_erase
        else 
            call psi_color (FA_COLOR(fa))
        
        # If the fill type is hollow, then just draw the polylines.
        # However, make sure the color of the polylines is the same as the
        # fill color and reset it after drawing.
        if (FA_STYLE(fa) == GF_HOLLOW || PSI_DOFILL(gkt) == NO) {
            pl = PSI_PLAP(g_kt)
            old_pl_ltype = PL_LTYPE(pl)
            if (FA_STYLE(fa) == GF_CLEAR)
                PL_LTYPE(pl) = GL_CLEAR
            old_pl_color = PL_COLOR(pl)
            PL_COLOR(pl) = FA_COLOR(fa)
            call psi_polyline (p, npts)
            PL_COLOR(pl) = old_pl_color
            PL_LTYPE(pl) = old_pl_ltype
        }
        
        # Else, set the fill style and path.
        else {
            
            # To isolate the effects of the fill area routine, suround it
            # with a gsave-grestore.  The grestore is done in the FI operator.
            call psk_out ("GS")
            
            # If the style is GF_CLEAR, GF_HOLLOW, or GF_SOLID, no pattern needs
            # to be defined.  However, any other pattern is specified by
            # HATCHn where n is the style index.
            if (FA_STYLE(fa) != GF_CLEAR &&
                FA_STYLE(fa) != GF_HOLLOW &&
                FA_STYLE(fa) != GF_SOLID) {
                call sprintf (g_output, SZ_LINE, "H%d SP")
                call pargi (FA_STYLE(fa))
                call psk_out (g_output)
            }
            
            # Set the path.
            x = p[1]
            y = p[2]
            call psi_move (x, y)
            
            for (i=3;  i <= len_p;  i=i+2) {
                x = p[i]
                y = p[i+1]
                call psi_draw (x, y)
            }
            
            # Close the path, fill, and restore the graphics context.
            call psk_out ("FI")
	    g_nseg = 0
            
        }
        
end
#---------------------------------------------------------------------------
# End of psi_fillarea
#---------------------------------------------------------------------------

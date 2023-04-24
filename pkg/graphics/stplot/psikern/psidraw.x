include "psi.h"

#---------------------------------------------------------------------------
.help psi_draw 1May92 plot
.ih
NAME
psi_draw -- Output a pen draw instruction.
.ih
USAGE
call psi_draw (x, y)
.ih
ARGUMENTS
.ls x y (int)
The NDC point to draw to.
.le
.ih
DESCRIPTION
The draw command is implemented in PostScript by the psikern-defined
command 'D'.
.ih
SEE ALSO
t_psikern, psk_prolog, psi_prolog.ps
.endhelp
#---------------------------------------------------------------------------

procedure psi_draw (x, y)

int     x, y    # I:  point to draw to

# Declarations.
int sx, sy      # Coordinates in output space.

pointer coords  # String representing the compressed coordinates.

include "psi.com"

# Function prototypes.
pointer psi_compress_coords()

begin
        # Set how many drawing operations have been done.
        g_ndraw = g_ndraw + 1

	# Set DCS page indicator.
	if (g_bpage) {
	    call psk_page
	    g_bpage = false
	}

        # Set how many line segments have been drawn.
        if (g_nseg > PS_MAX_SEGMENTS)
            call psi_stroke

        # Convert and compress the coordinates.
        sx = x * PSI_GKI2OUT(g_kt)
        sy = y * PSI_GKI2OUT(g_kt)
        coords = psi_compress_coords (sx, sy)
        
        # Output the instruction.
        call sprintf (g_output, SZ_LINE, "D %s")
        call pargstr (Memc[coords])
        call psk_out (g_output)
        g_nseg = g_nseg + 1
        
        # Free memory.
        call mfree (coords, TY_CHAR)
        
end
#---------------------------------------------------------------------------
# End of psi_draw
#---------------------------------------------------------------------------

include "psi.h"

#---------------------------------------------------------------------------
.help psi_move 1May92 plot
.ih
NAME
psi_move -- Output a pen move instruction.
.ih
USAGE
call psi_move (x, y)
.ih
ARGUMENTS
.ls x, y (int)
The NDC point to move to.
.le
.ih
DESCRIPTION
The move instruction is implemented by the PSI-defined PostScript command
'M'.
.ih
SEE ALSO
t_psikern, psk_prolog
.endhelp
#---------------------------------------------------------------------------

procedure psi_move (x, y)

int     x, y                    # point to move to

# Declarations.
int sx, sy       # The scaled coordinates.

pointer coords   # The string representation of the coordinates.

include "psi.com"

# Function prototypes.
pointer psi_compress_coords()

begin
	# Set DCS page indicator.
	if (g_bpage) {
	    call psk_page
	    g_bpage = false
	}
        
        # Convert and compress the coordinates.
        sx = x * PSI_GKI2OUT(g_kt)
        sy = y * PSI_GKI2OUT(g_kt)
        coords = psi_compress_coords (sx, sy)
        
        # Move the point in PostScript
        call sprintf (g_output, SZ_LINE, "M %s")
        call pargstr (Memc[coords])
        call psk_out (g_output)
        
        # Free the string.
        call mfree (coords, TY_CHAR)
        
end
#---------------------------------------------------------------------------
# End of psi_move
#---------------------------------------------------------------------------

include <gki.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_linewidth 1May92 plot
.ih
NAME
psi_linewidth -- Output a line width set instruction.
.ih
USAGE
call psi_linewidth (width)
.ih
ARGUMENTS
.ls width (real)
The new linewidth.
.le
.ih
DESCRIPTION
The line width is set by the PSI-defined PostScript command 'slw'.
.ih
SEE ALSO
t_psikern, psk_prolog
.endhelp
#---------------------------------------------------------------------------

procedure psi_linewidth (width)

int     width                   # new line width

# Declarations
real pswidth        # Line width in GKI/PostScript units.

include "psi.com"

begin
	# Set DCS page indicator.
	if (g_bpage) {
	    call psk_page
	    g_bpage = false
	}
        
        # See if the width has changed, and if so, tell the device.
        if (PSI_WIDTH(g_kt) != width) {
            pswidth = GKI_UNPACKREAL(width)
            if (pswidth >= 1)
                pswidth = PSI_LWORIGIN(g_kt) +
                           (PSI_LWSLOPE(g_kt) * PSI_LWORIGIN(g_kt) *
                             (pswidth - 1))
            else
                pswidth = PSI_LWORIGIN(g_kt) -
                           (PSI_LWSLOPE(g_kt) * PSI_LWORIGIN(g_kt) *
                            (1 - pswidth))
            
            call sprintf (g_output, SZ_LINE, "%d SW")
            call pargr (pswidth)
            call psk_out (g_output)
            
            PSI_WIDTH(g_kt) = width
        }
        
end
#---------------------------------------------------------------------------
# End of psi_linewidth
#---------------------------------------------------------------------------

include "psi.h"

#---------------------------------------------------------------------------
.help psi_erase 1May92 plot
.ih
NAME
psi_erase -- Set up so the following drawing "erases" the page.
.ih
USAGE
call psi_erase
.ih
DESCRIPTION
In PostScript, there is no concept of erasure.  However, a similar effect
can be obtained by painting in the appropriate "color", the background.
In general, PostScript devices have a white background- full color.
.endhelp
#---------------------------------------------------------------------------
procedure psi_erase()

# Declarations.
include "psi.com"

begin
        # Set DCS page indicator.
	if (g_bpage) {
	    call psk_page
	    g_bpage = false
	}
	
        # Set the appropriate color.
        call psk_out ("1 1 1 setrgbcolor")
        
        # Set that the current color is undefined.
        PSI_COLOR(g_kt) = -1
        
end
#---------------------------------------------------------------------------
# End of psi_erase
#---------------------------------------------------------------------------

include	"psk.h"

#---------------------------------------------------------------------------
.help psk_page Oct93 source
.ih
NAME
psk_page -- Indicate in the output that a new page has begun.
.endhelp
#---------------------------------------------------------------------------
procedure psk_page()

# Declarations
include	"psk.com"

begin
        # Output the document page for the next page.
        call psk_flush
        call sprintf (ps_output, MAX_CHAR, "%%%%Page: %d %d")
        call pargi (ps_frame + 1)
        call pargi (ps_frame + 1)
        call psk_out (ps_output)
        call psk_flush
	
	call sprintf (ps_output, MAX_CHAR, "NP")
	call psk_out (ps_output)
	call psk_flush
end
#---------------------------------------------------------------------------
# End of psk_page
#---------------------------------------------------------------------------

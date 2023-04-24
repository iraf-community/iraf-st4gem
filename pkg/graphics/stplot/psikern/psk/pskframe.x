include "psk.h"

#---------------------------------------------------------------------------
.help psk_frame 1May92 plot
.ih
NAME
psk_frame -- Output a page instruction.
.ih
USAGE
call psk_frame
.ih
ARGUMENTS
.ih
DESCRIPTION
For want of a better definition, a GKI "frame" is a PostScript "page"
Thus to change frames is to output the current page.  This definition may
need revision for Display PostScript.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psk_frame

# Declarations.
include "psk.com"

begin
        
        if (ps_debug)
            call eprintf ("psk_frame: start a new frame\n")
        
        # Output the page command.
        call psk_out ("CL")
        
        # Set indicators.
        ps_frame = ps_frame + 1
end
#---------------------------------------------------------------------------
# End of psk_frame
#---------------------------------------------------------------------------

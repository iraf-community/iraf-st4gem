include "psi.h"

#---------------------------------------------------------------------------
.help psi_stroke 1May92 plot
.ih
NAME
psi_stroke -- Make sure the last set of drawing commands are set on paper.
.ih
USAGE
call psi_stroke
.ih
DESCRIPTION
In PostScript, drawing operations are not "set" until the current path
is finished.  Finish the current path and place it on the page.
.ih
SEE ALSO
t_psikern, psk_prolog
.endhelp
#---------------------------------------------------------------------------

procedure psi_stroke

# Declarations.
include "psi.com"

begin
        
        call psk_out ("SR")
        g_nseg = 0
        
end
#---------------------------------------------------------------------------
# End of psi_stroke
#---------------------------------------------------------------------------

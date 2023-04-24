include "psi.h"

#---------------------------------------------------------------------------
.help psi_cancel 1May92 plot
.ih
NAME
psi_cancel -- Cancel any buffered output.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_cancel (dummy)

int     dummy                   # not used at present

# Declarations.
include "psi.com"

begin
        
        if (g_debug)
            call eprintf ("psi_cancel: Cancel current operation.\n")
        
        if (g_kt == NULL)
            return
        call psi_reset
        
end
#---------------------------------------------------------------------------
# End of psi_cancel
#---------------------------------------------------------------------------

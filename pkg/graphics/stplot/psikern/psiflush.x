include "psi.h"

#---------------------------------------------------------------------------
.help psi_flush 1May92 plot
.ih
NAME
psi_flush -- Flush the output file.
.ih
USAGE
call psi_flush (dummy)
.ih
ARGUMENTS
.ls dummy (int)
Not used, only needed for the GKI interface.
.le
.ih
DESCRIPTION
Doesn't do much but call psk_flush (which doesn't do much either).
.ih
SEE ALSO
t_psikern, psk_flush
.endhelp
#---------------------------------------------------------------------------

procedure psi_flush (dummy)

int dummy  # Place holder, nothing more

# Declarations.
include "psi.com"

begin
        
        if (g_debug)
            call eprintf ("psi_flush: Flushing current buffered output.\n")
        
        call psk_flush
        
end
#---------------------------------------------------------------------------
# End of psi_flush
#---------------------------------------------------------------------------

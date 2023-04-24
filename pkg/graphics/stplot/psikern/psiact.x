include "psi.h"

#---------------------------------------------------------------------------
.help psi_[de,re]actws Dec92 source
.ih
NAME
psi_deactws, psi_reactws -- Deactivate and reactivate the output device.
.endhelp
#---------------------------------------------------------------------------
procedure psi_deactws

# Declarations.
include "psi.com"

begin
        if (g_debug)
            call eprintf ("psi_deactws: Deactivating\n")

        call psk_flush
        call psk_out ("%PSKCloseWorkStation")
        call psk_flush
end
#---------------------------------------------------------------------------
# End of psi_deactws
#---------------------------------------------------------------------------
procedure psi_reactws

# Declarations.
include "psi.com"

begin
        if (g_debug)
            call eprintf ("psi_reactws: Reactivating\n")

        call psk_flush
        call psk_out ("%PSKOpenWorkStation")
        call psk_flush
end
#---------------------------------------------------------------------------
# End of psi_reactws
#---------------------------------------------------------------------------

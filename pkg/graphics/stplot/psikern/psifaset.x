include <gki.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_faset 1May92 plot
.ih
NAME
psi_faset -- Set the fillarea attributes.
.ih
USAGE
call psi_faset (gki)
.ih
ARGUMENTS
.ls gki (short[ARB])
The GKI attribute structure.
.le
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_faset (gki)

short   gki[ARB]                # attribute structure

# Declarations.
pointer fa

include "psiparams.com"
include "psi.com"

begin
        
        if (g_debug)
            call eprintf ("psi_faset: Setting fillarea attributes.\n")
        
        fa = PSI_FAAP(g_kt)
        FA_STYLE(fa) = gki[GKI_FASET_FS]
        FA_COLOR(fa) = gki[GKI_FASET_CI]
        if (!IS_INDEFI(area_color) && FA_COLOR(fa) == DEF_COLOR)
            FA_COLOR(fa) = area_color
        
end
#---------------------------------------------------------------------------
# End of psi_faset
#---------------------------------------------------------------------------

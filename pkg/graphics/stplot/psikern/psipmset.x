include <gki.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_pmset 1May92 plot
.ih
NAME
psi_pmset -- Set the polymarker attributes.
.ih
USAGE
call psi_pmset (gki)
.ih
ARGUMENTS
.ls gki (short[ARB])
The GKI attribute structure
.le
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_pmset (gki)

short   gki[ARB]                # attribute structure

# Declarations
pointer pm

include "psiparams.com"
include "psi.com"

begin
        
        if (g_debug)
            call eprintf ("psi_pmset: Setting polymarker attributes.\n")
        
        pm = PSI_PMAP(g_kt)
        PM_LTYPE(pm) = gki[GKI_PMSET_MT]
        PM_WIDTH(pm) = gki[GKI_PMSET_MW]
        PM_COLOR(pm) = gki[GKI_PMSET_CI]
        if (!IS_INDEFI(marker_color) && PM_COLOR(pm) == DEF_COLOR)
            PM_COLOR(pm) = marker_color
        
end
#---------------------------------------------------------------------------
# End of psi_pmset
#---------------------------------------------------------------------------

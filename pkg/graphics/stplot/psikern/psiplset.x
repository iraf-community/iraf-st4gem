include <gki.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_plset 1May92 plot
.ih
NAME
psi_plset -- Set the polyline attributes.  
.ih
USAGE
call psi_plset (gki)
.ih
ARGUMENTS
.ls gki (short[ARB])
The attribute structure.
.le
.ih
DESCRIPTION
The polyline width parameter is passed to the encoder as a packed floating 
point number, i.e., int(LWx100).
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_plset (gki)

short   gki[ARB]                # attribute structure
pointer pl

# Declarations.
include "psiparams.com"
include "psi.com"

begin
        
        if (g_debug)
            call eprintf ("psi_plset: Setting polyline attributes.\n")
        
        pl = PSI_PLAP(g_kt)
        PL_LTYPE(pl) = gki[GKI_PLSET_LT]
        PL_WIDTH(pl) = gki[GKI_PLSET_LW]
        PL_COLOR(pl) = gki[GKI_PLSET_CI]
        if (!IS_INDEFI(line_color) && PL_COLOR(pl) == DEF_COLOR)
            PL_COLOR(pl) = line_color
        
end
#---------------------------------------------------------------------------
# End of psi_plset
#---------------------------------------------------------------------------

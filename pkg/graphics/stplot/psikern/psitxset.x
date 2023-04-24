include <gset.h>
include <gki.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_txset 1May92 plot
.ih
NAME
psi_txset -- Set the text drawing attributes.
.ih
USAGE
call psi_txset (gki)
.ih
ARGUMENTS
.ls gki (short[ARB])
GKI Attribute structure.
.le
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_txset (gki)

short   gki[ARB]                # attribute structure

# Declarations.
pointer tx

include "psiparams.com"
include "psi.com"

begin
        
        if (g_debug)
            call eprintf ("psi_txset: Setting text attributes.\n")
        
        tx = PSI_TXAP(g_kt)
        
        TX_UP(tx)       = gki[GKI_TXSET_UP] 
        TX_PATH(tx)     = gki[GKI_TXSET_P ] 
        TX_HJUSTIFY(tx) = gki[GKI_TXSET_HJ] 
        TX_VJUSTIFY(tx) = gki[GKI_TXSET_VJ] 
        TX_FONT(tx)     = gki[GKI_TXSET_F ]
        TX_QUALITY(tx)  = gki[GKI_TXSET_Q ] 
        TX_COLOR(tx)    = gki[GKI_TXSET_CI]
        if (!IS_INDEFI(text_color) && TX_COLOR(tx) == DEF_COLOR)
            TX_COLOR(tx) = text_color
        
        TX_SPACING(tx)  = GKI_UNPACKREAL(gki[GKI_TXSET_SP])
        TX_SIZE(tx)     = gki[GKI_TXSET_SZ]
        
end
#---------------------------------------------------------------------------
# End of psi_txset
#---------------------------------------------------------------------------

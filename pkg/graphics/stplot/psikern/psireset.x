include <gki.h>
include <gset.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_reset 1May92 plot
.ih
NAME
PSI_RESET -- Reset the state of the transform common.
.ih
USAGE
call psi_reset()
.ih
DESCRIPTION
Reset the state of the transform common, i.e., in response to
a clear or a cancel.  Initialize all attribute packets to their default
values and set the current state of the device to undefined, forcing the
device state to be reset when the next output instruction is executed.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_reset

# Declarations.
pointer pl, pm, fa, tx

include "psi.com"

begin
	if (g_debug)
	    call eprintf ("psi_reset: Resetting internal parameters.\n")
	
        # Set the maximum color back to the original.
        g_maxgrcolor = g_defmaxgrcolor
        
        # Set pointers to attribute substructures.
        pl = PSI_PLAP(g_kt)
        pm = PSI_PMAP(g_kt)
        fa = PSI_FAAP(g_kt)
        tx = PSI_TXAP(g_kt)
        
        # Initialize the attribute packets.
        PL_LTYPE(pl)    = 1
        PL_WIDTH(pl)    = GKI_PACKREAL(1.)
        PL_COLOR(pl)    = 1
        PM_LTYPE(pm)    = 1
        PM_WIDTH(pm)    = GKI_PACKREAL(1.)
        PM_COLOR(pm)    = 1
        FA_STYLE(fa)    = 1
        FA_COLOR(fa)    = 1
        TX_UP(tx)       = 90
        TX_SIZE(tx)     = GKI_PACKREAL(1.)
        TX_PATH(tx)     = GT_RIGHT
        TX_HJUSTIFY(tx) = GT_LEFT
        TX_VJUSTIFY(tx) = GT_BOTTOM
        TX_FONT(tx)     = GT_ROMAN
        TX_COLOR(tx)    = 1
        TX_SPACING(tx)  = 0.0
        TX_VARIABLE(tx) = PSI_DEFVAR(g_kt)
        
        # Set the device attributes to undefined, forcing them to be reset
        # when the next output instruction is executed.
        PSI_TYPE(g_kt)          = -1
        PSI_WIDTH(g_kt)         = -1
        PSI_COLOR(g_kt)         = -1
        PSI_TXSIZE(g_kt)        = -1
        PSI_TXFONT(g_kt)        = -1
        PSI_UP(g_kt)            = INDEFI
        PSI_VARIABLE(g_kt)      = -1
        PSI_UP(g_kt)            = INDEFI
        PSI_PATH(g_kt)          = -1
        PSI_HJUSTIFY(g_kt)      = -1
        PSI_VJUSTIFY(g_kt)      = -1
        PSI_SPACING(g_kt)       = INDEFR
        
end
#---------------------------------------------------------------------------
# End of psi_reset
#---------------------------------------------------------------------------

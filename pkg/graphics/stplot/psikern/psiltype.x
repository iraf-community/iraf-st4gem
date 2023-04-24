include <gset.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_ltype 1May92 plot
.ih
NAME
psi_ltype - Set the polyline type in PostScript.
.ih
USAGE
call psi_ltype (ltype)
.ih
ARGUMENTS
.ls ltype (int)
I: The line type.  If the line type hasn't changed, don't output any commands
to change it.
.le
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_ltype (ltype)

int ltype  # I:  The line type to set.

# Declarations
include "psi.com"

begin
        # Set DCS page indicator.
	if (g_bpage) {
	    call psk_page
	    g_bpage = false
	}

	# Set the line type.
        if (PSI_TYPE(g_kt) != ltype)
            switch (ltype) {
            case GL_SOLID:
                call psk_out ("LPS LT")
                
            case GL_DASHED:
                call psk_out ("LPD LT")
                
            case GL_DOTTED:
                call psk_out ("LPP LT")
                
            case GL_DOTDASH:
                call psk_out ("LDD LT")
            }
        
        PSI_TYPE(g_kt) = ltype
        
end
#---------------------------------------------------------------------------
# End of psi_ltype
#---------------------------------------------------------------------------

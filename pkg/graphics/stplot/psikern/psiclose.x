include "psi.h"

#---------------------------------------------------------------------------
.help psi_close 1May92 st4gem.stplot
.ih
NAME
psi_close -- Close the PSI translation kernel.  
.ih
USAGE
call psi_close
.ih
DESCRIPTION
Close the PSI translation kernel.  Close the spool file so
the output is finally plotted.  Free up storage.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_close()

# Declarations.
include "psi.com"

begin
        
        if (g_debug)
            call eprintf ("psi_close: Closing everything down.\n")
        
        # Close output metafile, disposing of it to the host system.
        call psi_reactws
        call psk_close ((g_ndraw > 0 || g_nframes > 0),
			PSI_FONTLIST(g_kt))
        
        # Close the prolog file.
        call close (PSI_PROLOG(g_kt))
        
        # Free kernel data structures.
        call mfree (PSI_SBUF(g_kt), TY_CHAR)
        call mfree (g_kt, TY_STRUCT)
        g_kt = NULL
        
end
#---------------------------------------------------------------------------
# End of psi_close
#---------------------------------------------------------------------------

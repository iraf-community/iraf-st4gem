include "psi.h"

#---------------------------------------------------------------------------
.help psi_closews 1May92 plot
.ih
NAME
PSI_CLOSEWS -- Close the named workstation.  
.ih
USAGE
call psi_closews (devname, n)
.ih
ARGUMENTS
.ls devname (short[n])
Device name (not used)
.le
.ls n (int)
Number of characters in device name.
.le
.ih
DESCRIPTION
This is a noop.  There is nothing to flush, etc. at the moment.
The spool file is closed only on the next plot or at gktclose time.
If the spool file is closed here, APPEND mode would not work.
.ih
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_closews (devname, n)

short   devname[ARB]            # device name 
int     n                       # length of device name

pointer sp                      # Stack pointer.
pointer xstr                    # Generic.

# Declarations.
include "psi.com"

begin
        if (g_debug) {
            call smark (sp)
            call salloc (xstr, n, TY_CHAR)

            call achtsc (devname, Memc[xstr], n)
            Memc[xstr+n] = EOS
            call sprintf (g_output, SZ_LINE,
			  "psi_closews: Closing device %s.\n")
            call pargstr (Memc[xstr])
            call eprintf (g_output)
            
            call sfree (sp)
        }
	if (g_maxframes < 0)
	    call psi_clear (0)

        call psk_flush
        call psk_out ("%PSKCloseWorkStation")
        call psk_flush
        
end
#---------------------------------------------------------------------------
# End of psi_closews
#---------------------------------------------------------------------------

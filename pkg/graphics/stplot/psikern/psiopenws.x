include <error.h>
include <gki.h>
include <mach.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_openws 1May92 plot
.ih
NAME
psi_openws -- Open the named workstation.  
.ih
USAGE
call psi_openws (devname, n, mode)
.ih
ARGUMENTS
.ls devname (short[ARB])
Device name packed into short integers.
.le
.ls n (int)
Number of characters in devname
.le
.ls mode (int)
Access mode to open the devname file.
.le
.ih
DESCRIPTION
Once a workstation has been opened we leave it open until some other 
workstation is opened or the kernel is closed.  Opening a workstation 
involves initialization of the kernel data structures, following by 
initialization of the device itself.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_openws (devname, n, mode)

short   devname[ARB]            # device name
int     n                       # length of device name
int     mode                    # access mode

# Declarations.
bool need_open  # True if the device needs to be opened.
bool same_dev   # True if this device was used in the last GKI operations.

pointer buf     # Temporary buffer.
pointer sp      # Stack pointer

include "psi.com"
include "psiparams.com"

# Function prototypes.
pointer ttygdes ()
bool    streq ()

begin
        
        if (g_debug)
            call eprintf ("psi_openws: Opening the workstation.\n")
        
        # Allocate memory.
        call smark (sp)
        call salloc (buf, max (SZ_FNAME, n), TY_CHAR)
        
        # If a device was named when the kernel was opened then output will
        # always go to that device (g_device) regardless of the device named
        # in the OPENWS instruction.  If no device was named (null string)
        # then unpack the device name, passed as a short integer array.
        if (g_device[1] == EOS) {
            call achtsc (devname, Memc[buf], n)
            Memc[buf+n] = EOS
        } else
            call strcpy (g_device, Memc[buf], SZ_FNAME)
 
#	call eprintf("devname was %d, while device is set to %s\n")
#		call pargi(devname)
#		call pargstr(g_device)


        # Find out if first time, and if not, if same device as before
        # note that if (g_kt == NULL), then same_dev is false.
        same_dev = false
        need_open = true
        
        if (g_kt != NULL) {
            same_dev = (streq (Memc[PSI_DEVNAME(g_kt)], Memc[buf]))
            if (!same_dev)
                call psk_close ((g_ndraw > 0 || g_nframes > 0),
				PSI_FONTLIST(g_kt))
            else
                need_open = false
        }
        
        # Initialize the kernel data structures.  Open graphcap descriptor
        # for the named device, allocate and initialize descriptor and common.
        # graphcap entry for device must exist.
        if (need_open) {
            if (!same_dev) {
                if (g_kt != NULL)
                    call ttycdes (g_tty)
                iferr (g_tty = ttygdes (Memc[buf]))
                    call erract (EA_ERROR)
                

#        call eprintf("Buf = %s, tty = %d\n")
#                call pargstr(Memc[buf])
#		call pargi(g_tty)

                # Initialize data structures if we had to open a new device.
                call psi_init (g_tty, Memc[buf])
                call psi_reset
            }
            
            # Open the output file.  Metacode output to the device will be
            # spooled and then disposed of to the device at CLOSEWS time.
            iferr (call psk_open (Memc[PSI_DEVNAME(g_kt)], g_tty, out_file)) {
                call ttycdes (g_tty)
                call erract (EA_ERROR)
            } else {
                
                # Send the prolog.
                call psi_prolog
                
                # Reset the statistics.
                g_nframes = 0
                g_ndraw = 0
		g_bpage = true
            }
        }

        # Else, just send an indicator that the workstation is reopening.
        else {
            call psk_flush
            call psk_out ("%PSKOpenWorkStation")
            call psk_flush
        }
        
        # Clear the screen if device is being opened in new_file mode.
        # This is a nop if we really opened a new device, but it will clear
        # the screen if this is just a reopen of the same device in new file
        # mode.
        if (mode == NEW_FILE)
            call psi_clear (0)
        
        # Release memory.
        call sfree (sp)
        
end
#---------------------------------------------------------------------------
# End of psi_openws
#---------------------------------------------------------------------------

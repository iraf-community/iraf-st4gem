include <mach.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_clear 1May92 stsdast.stplot
.ih
NAME
PSI_CLEAR -- Start a new page.
.ih
USAGE
call psi_clear (dummy)
.ih
ARGUMENTS
.ls dummy (int)
Not used at present.
.le
.ih
DESCRIPTION
Start a new PostScript page..  All attribute packets are
initialized to their default values.  Redundant calls or calls immediately
after a workstation open (before anything has been drawn) are ignored.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_clear (dummy)

int     dummy                   # not used at present

# Declarations.
include "psi.com"

errchk  psk_open

begin
        if (g_debug)
            call eprintf ("psi_clear: Clearing device.\n")

        # This is a no-op if nothing has been drawn.
        if (g_kt == NULL || g_ndraw == 0)
            return
        
        # Start a new frame.  This is done either by issuing the frame advance
        # instruction or by starting a new metafile.  Close the output file and 
        # start a new metafile if the maximum frame count has been reached.
        # This disposes of the metafile to the system, causing the actual
        # plots to be drawn.  Open a new metafile ready to receive next frame.
        g_nframes = g_nframes + 1
        if (g_nframes >= g_maxframes) {
            call psk_close (true, PSI_FONTLIST(g_kt))
            call psk_open (Memc[PSI_DEVNAME(g_kt)], g_tty, "")
            call psi_prolog
            g_nframes = 0
        } else
            call psk_frame
        
        # Init kernel data structures.
        call psi_reset
        g_ndraw = 0
	g_bpage = true
end
#---------------------------------------------------------------------------
# End of psi_clear
#---------------------------------------------------------------------------

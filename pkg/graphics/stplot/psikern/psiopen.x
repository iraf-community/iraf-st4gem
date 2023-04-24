include <gki.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_open 1May92 plot
.ih
NAME
psi_open -- Install the PSI kernel as a graphics kernel device driver.
.ih
USAGE
call psi_open (devname, dd)
.ih
ARGUMENTS
.ls devname (char[ARB])
Nonnull for forced output to a device.
.le
.ls dd (int[ARB])
Device table to be initialized.
.le
.ih
DESCRIPTION
The device table DD consists of an array of the entry point addresses for
the driver procedures.  If a driver does not implement a particular
instruction the table entry for that procedure may be set to zero, causing
the interpreter to ignore the instruction.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_open (devname, dd)

char    devname[ARB]            # nonnull for forced output to a device
int     dd[ARB]                 # device table to be initialized

# Declarations.
int len_devname  # Lenght of the device name.

pointer devns    # Device name.
pointer sp       # Stack pointer.

include "psi.com"

# Function prototypes.
int locpr(), strlen()

extern psi_cancel(), psi_clear(), psi_closews(), psi_deactws()
extern psi_escape(), psi_faset(), psi_fillarea(), psi_flush()
extern psi_openws(), psi_plset(), psi_pmset(), psi_polyline()
extern psi_polymarker(), psi_putcellarray(), psi_reactws()
extern psi_text(), psi_txset()

begin
        
        # Allocate memory.
        call smark (sp)
        call salloc (devns, SZ_FNAME, TY_SHORT)
        
        # Flag first pass.  Save forced device name in common for OPENWS.
        # Zero the frame and instruction counters.
        g_kt = NULL
        call strcpy (devname, g_device, SZ_GDEVICE)
        
        # Install the device driver.
        dd[GKI_OPENWS]          = locpr (psi_openws)
        dd[GKI_CLOSEWS]         = locpr (psi_closews)
        dd[GKI_DEACTIVATEWS]    = locpr (psi_deactws)
        dd[GKI_REACTIVATEWS]    = locpr (psi_reactws)
        dd[GKI_MFTITLE]         = 0                     # Not Applicable.
        dd[GKI_CLEAR]           = locpr (psi_clear)
        dd[GKI_CANCEL]          = locpr (psi_cancel)
        dd[GKI_FLUSH]           = locpr (psi_flush)
        dd[GKI_POLYLINE]        = locpr (psi_polyline)
        dd[GKI_POLYMARKER]      = locpr (psi_polymarker)
        dd[GKI_TEXT]            = locpr (psi_text)
        dd[GKI_FILLAREA]        = locpr (psi_fillarea)
        dd[GKI_PUTCELLARRAY]    = locpr (psi_putcellarray)
       # dd[GKI_PUTCELLARRAY]    = 0
        dd[GKI_SETCURSOR]       = 0                     # Not Applicable.
        dd[GKI_PLSET]           = locpr (psi_plset)
        dd[GKI_PMSET]           = locpr (psi_pmset)
        dd[GKI_TXSET]           = locpr (psi_txset)
        dd[GKI_FASET]           = locpr (psi_faset)
        dd[GKI_GETCURSOR]       = 0                     # Not Applicable.
        dd[GKI_GETCELLARRAY]    = 0                     # Not Applicable.
        dd[GKI_ESCAPE]          = locpr (psi_escape)
        dd[GKI_SETWCS]          = 0                     # Not Applicable.
        dd[GKI_GETWCS]          = 0                     # Not Applicable.
        dd[GKI_UNKNOWN]         = 0                     # Proper action
        
        # If a device was named open the workstation as well.  This is
        # necessary to permit processing of metacode files which do not
        # contain the open workstation instruction.
        len_devname = strlen (devname)
        if (len_devname > 0) {
            call achtcs (devname, Mems[devns], len_devname)
            call psi_openws (Mems[devns], len_devname, NEW_FILE)
        }

#	call eprintf("PSIOPEN: device name %s packed into %d\n")
#		call pargstr(devname)
#		call pargi(Mems[devns])
        
        # Release memory.
        call sfree (sp)
        
end
#---------------------------------------------------------------------------
# End of psi_open
#---------------------------------------------------------------------------

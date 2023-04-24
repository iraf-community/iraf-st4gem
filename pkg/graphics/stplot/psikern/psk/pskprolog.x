include <error.h>
include "psk.h"

#---------------------------------------------------------------------------
.help psk_prolog 1May92 plot
.ih
NAME
psk_prolog - Output the PostScript prolog
.ih
USAGE
call psk_prolog (prolog_fd)
.ih
ARGUMENTS
.ls prolog_fd (int)
The file descriptor of the PSIKERN PostScript prolog file.
.le
DESCRIPTION
Output to the spool file the PostScript prolog that contains the dictionary,
font, and variable setups necessary by the rest of the page description that
will be generated.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psk_prolog (prolog_fd)

int prolog_fd   # I:  The name of the PSIKERN prolog file.

# Declarations
include "psk.com"

# Function prototypes.
int errcode()

begin
        
        if (ps_debug)
            call eprintf ("psk_prolog: outputting prolog.\n")
        
        # Flush what buffered output there might be.
        call psk_flush
        
        # Copy the prolog file onto the output PostScript file.
        call fprintf (ps_fd, "\n")
        iferr (call fcopyo (prolog_fd, ps_fd)) {
            call eprintf ("psk_prolog: Could not copy prolog onto spool file.\n")
            call erract (errcode())
        }
        
        # Set back to the beginning of the prolog file for the event that
        # it will need to be copied again.
        call seek (prolog_fd, BOF)
        
end
#---------------------------------------------------------------------------
# End of psk_prolog
#---------------------------------------------------------------------------

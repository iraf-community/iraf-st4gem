include "psk.h"

#---------------------------------------------------------------------------
.help psk_flush 1May92 plot
.ih
NAME
psk_flush -- Flush the output to a file.
.ih
USAGE
call psk_flush
.ih
ARGUMENTS
.ih
DESCRIPTION
At the moment, this is just a call to flush().
.ih
SEE ALSO
t_psikern, psi_polyline
.endhelp
#---------------------------------------------------------------------------

procedure psk_flush

# Declarations
include "psk.com"

begin
        
        # Dump out the rest of the output buffer.
        if (ps_buffer_eos > 1 && ps_buffer_output){
            ps_buffer[ps_buffer_eos] = EOS
            call fprintf (ps_fd, "%s\n")
            call pargstr (ps_buffer)
            ps_buffer_eos = 1
        }
        
        # System-wise flush the stream.
        call flush (ps_fd)
        
end
#---------------------------------------------------------------------------
# End of psk_flush
#---------------------------------------------------------------------------

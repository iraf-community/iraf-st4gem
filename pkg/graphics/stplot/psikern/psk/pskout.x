include "psk.h"

#---------------------------------------------------------------------------
.help psk_out 1May92 plot
.ih
NAME
psk_out -- Output the string to a buffer and flush the buffer if necessary.
.ih
USAGE
call psk_out (input)
.ih
ARGUMENTS
.ls input (char[ARB])
The string to buffer out.
.le
.ih
DESCRIPTION
This routine buffers output to the PostScript file to create a more compact
PostScript coding.  The goal is to minimize the number of newline characters
in the output.  A space is inserted after each input string in the buffer
to make sure commands aren't run together.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psk_out (input)

char input[ARB]         # I:  The string to buffer out.

# Declarations.
int input_length        # The number of characters in the input string.

include "psk.com"

# Function prototypes
int strlen()

begin
        
        # If not buffering, just write out the string.  Else, buffer.
        if (!ps_buffer_output) {
            call fprintf (ps_fd, "%s\n")
            call pargstr (input)
            
        } else {
            
            # If the buffer is full, flush it out.
            input_length = strlen (input)
            if ((input_length + 1 + ps_buffer_eos > MAX_CHAR) &&
                (ps_buffer_eos > 1)) {
                ps_buffer[ps_buffer_eos] = EOS
                call fprintf (ps_fd, "%s\n")
                call pargstr (ps_buffer)
                ps_buffer_eos = 1
            }
            
            # Copy the string into the buffer.  Check to make sure the input
            # string fits into the buffer.  If not, just write it out.
            if (input_length + 1 > MAX_CHAR) {
                call fprintf (ps_fd, "%s\n")
                call pargstr (input)
            } else {
                call strcpy (input, ps_buffer[ps_buffer_eos], MAX_CHAR)
                ps_buffer_eos = ps_buffer_eos + input_length
                ps_buffer[ps_buffer_eos] = ' '
                ps_buffer_eos = ps_buffer_eos + 1
                ps_buffer[ps_buffer_eos] = EOS
            }
            
        }
        
end
#---------------------------------------------------------------------------
# End of psk_out
#---------------------------------------------------------------------------

include <ctype.h>

#---------------------------------------------------------------------------
.help psi_mkpsstring 1May92 plot
.ih
NAME
psi_mkpsstring - Convert a string to a PostScript string.
.ih
USAGE
call psi_mkpsstring (input, output)
.ih
ARGUMENTS
.ls input (char[ARB])
The input string
.le
.ls output (char[ARB])
The PostScript string
.le
.ih
DESCRIPTION
This formats a string as a PostScript string literal.  The '(' and ')'
are prepended/appended and any embedded '(' and ')' are escaped.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_mkpsstring (input, output, maxchar)

char    input[ARB]   # I:  The string to convert.
char    output[ARB]  # O:  The PostScript string literal.
int     maxchar      # I:  The maximum number of characters in the output string.

# Declarations
int     ip              # Pointer to input array.
int     octal_ptr       # Pointer to octal char array.
int     op              # Pointer to output array.

char    octal[SZ_LINE]  # For conversion to octal representation.

begin
        
        output[1] = '('
        op = 2
        ip = 1
        while (input[ip] != EOS && op < maxchar - 2) {
            if (IS_PRINT(input[ip]))
                switch (input[ip]) {
                case '(', ')','\\':
                    output[op] = '\\'
                    output[op+1] = input[ip]
                    op = op + 2
                    
                default:
                    output[op] = input[ip]
                    op = op + 1
                }
            
            # Non-printable - format as octal code.
            else {
                call sprintf (octal, SZ_LINE, "%03o")
                call pargc (input[ip])
                octal_ptr = 1
                while (octal[octal_ptr] != EOS && op < maxchar -1) {
                    output[op] = octal[octal_ptr]
                    op = op +1
                    octal_ptr = octal_ptr + 1
                }
            }
            
            ip = ip + 1
        }
        
        output[op] = ')'
        output[op+1] = EOS
        
end
#---------------------------------------------------------------------------
# End of psi_mkpsstring
#---------------------------------------------------------------------------

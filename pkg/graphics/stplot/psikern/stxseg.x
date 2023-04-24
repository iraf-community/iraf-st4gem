include <gset.h>
include	"psi.h"

#---------------------------------------------------------------------------
.help stx_segment 1May92 plot
.ih
NAME
stx_segment -- Process the text string into segments.
.ih
USAGE
integer = stx_segment (text, n, out, start_font)
.ih
ARGUMENTS
.ls text (short[ARB])
The input packed text to decode.
.le
.ls n (int)
Number of characters in text.
.le
fs out (char[ARB])
The segmented output text string.  This has the form:
.nf
    out[1]           -- Font of the following segment.
    out[2]...out[n]  -- The text segment. 
    out[n+1]         -- EOS
    out[n+2]         -- Font for the next text segment.
     .....
    out[totlen - 1]  -- EOS
    out[totlen]      -- EOS and end of array of out.
.fi
.le
.ih
RETURNS
The total number of characters in out.
.ih
DESCRIPTION
Process the text string into segments, in the process
converting from type short to char.  The only text attribute that can
change within a string is the font, so segments are broken by font contols
etc. font select sequences embedded in the text.  The segments are encoded
sequentially in the output string.  The first character of each segment is
the font number.  A segment is delimited by EOS.  A font number of EOS
marks the end of the segment list.  The output string is assumed to be
large enough to hold the segmented text string.
.ih
HISTORY
.nf
   4Sep91 - Extracted from sys$gio/sgikern/sgitx.x and placed into its
            own file.  No modifications to the code, i.e. this is NOAO code.
            Jonathan D. Eisenhamer, STScI.
.fi
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

int procedure stx_segment (text, n, out, start_font)

short   text[ARB]               # input text
int     n                       # number of characters in text
char    out[ARB]                # output string
int     start_font              # initial font code

int     ip, op
int     totlen, font
          
begin
        out[1] = start_font
        totlen = 0
        op = 2
        
        for (ip=1;  ip <= n;  ip=ip+1) {
            if (text[ip] == '\\' && text[ip+1] == 'f') {
                # Select font.
                out[op] = EOS
                op = op + 1
                ip = ip + 2
                
                switch (text[ip]) {
                case 'B':
                    font = GT_BOLD
                case 'I':
                    font = GT_ITALIC
                case 'G':
                    font = GT_GREEK
		case 'P':
		    font = PS_USER
		case 'U':
		    font = PS_SETUP
                default:
                    font = GT_ROMAN
                }
                
                out[op] = font
                op = op + 1
                
            } else {
                # Deposit character in segment.
                out[op] = text[ip]
                op = op + 1
                totlen = totlen + 1
            }
        }
        
        # Terminate last segment and add null segment.
        
        out[op] = EOS
        out[op+1] = EOS
        
        return (totlen)
end
#---------------------------------------------------------------------------
# End of stx_segment
#---------------------------------------------------------------------------

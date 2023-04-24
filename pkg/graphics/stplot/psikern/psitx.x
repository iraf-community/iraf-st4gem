include <gset.h>
include <gki.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_text 1May92 plot
.ih
NAME
psi_text -- Draw a text string.  
.ih
USAGE
call psi_text (xc, yc, text, n)
.ih
ARGUMENTS
.ls xc, yc (int)
The point to draw the text string.
.le
.ls text (short[ARB])
The text string to draw.
.le
.ls n (int)
The number of characters in the text string.
.le
.ih
DESCRIPTION
The string is drawn at the position (X,Y) using the text attributes
set by the last GKI_TXSET instruction.  The text string to be drawn
may contain embedded set font escape sequences of the form \fR (roman), 
\fG (greek), etc.  We break the input text sequence up into segments at
font boundaries and draw these on the output device, setting the text
size, color, font, and position at the beginning of each segment.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_text (xc, yc, text, n)

int     xc, yc                  # I:  Where to draw text string
short   text[ARB]               # I:  Text string
int     n                       # I:  Number of characters

# Declarations.
real    ch, cw          # Character height, width.  
real    scale           # Scale size of characters.

int     font            # The current font.
int     nchar           # Total number of characters in segmented string (used
                        # only for return from stx_segment function)

pointer ip              # Index pointer into segmented string.
pointer out_text        # The text segment to be output.
pointer out_text2       # The text segment to be output.
pointer ps_text         # PostScript version of the text line.
pointer sp              # Stack pointer.
pointer tx              # GKI text attribute structure.

include "psi.com"

# Function prototypes.
int     stx_segment()
bool    itob()

begin
        
        if (g_debug){
            call eprintf ("psi_text: Writing text of %d characters.\n")
            call pargi (n)
        }
        
        call smark (sp)
        
        # Keep track of the number of drawing instructions since the last frame
        # clear.
        g_ndraw = g_ndraw + 1

	# Set DCS page indicator.
	if (g_bpage) {
	    call psk_page
	    g_bpage = false
	}
        
        # Set pointer to the text attribute structure.
        tx = PSI_TXAP(g_kt)
        
        # If the UP vector has changed, set it.
        if (PSI_UP(g_kt) != TX_UP(tx)) {
            call sprintf (g_output, SZ_LINE, "%d PA")
            call pargi (TX_UP(tx) - 90)
            call psk_out (g_output)
            PSI_UP(g_kt) = TX_UP(tx)
        }
        
        # Set text size if changed.
        if (PSI_TXSIZE(g_kt) != TX_SIZE(tx)) {
            scale = GKI_UNPACKREAL(TX_SIZE(tx))

            # Because fonts in PostScript tend to get small quickly,
            # adjust by using the root of the size.  The same problem
            # holds for large fonts, but it is not nearly so noticable.
            # The crossing point, .5, is arbitrary and looked like a reasonable
            # place to go from linear to root.  The larger this point is
            # the larger small fonts will be.  At 1 however, small fonts
            # were a bit too large.
            if (scale < .5)
                scale = 0.7071068 * sqrt (scale)
            ch = PSI_CHARHEIGHT(g_kt,1) * scale * PSI_GKI2OUT(g_kt) 
            cw = PSI_CHARWIDTH(g_kt,1) * scale * PSI_GKI2OUT(g_kt) 
            call sprintf (g_output, SZ_LINE, "%g %g FS")
            call pargr (cw)
            call pargr (ch)
            call psk_out (g_output)
            PSI_TXSIZE(g_kt) = TX_SIZE(tx)
        }
        
        # Set path if changed.
        if (PSI_PATH(g_kt) != TX_PATH(tx)) {
            if (TX_PATH(tx) == GT_RIGHT || TX_PATH(tx) == GT_LEFT) {
                call sprintf (g_output, SZ_LINE, "%d PH")
                call pargi (GT_RIGHT)
            } else {
                call sprintf (g_output, SZ_LINE, "%d PH")
                call pargi (GT_UP)
            }
            call psk_out (g_output)
            PSI_PATH(g_kt) = TX_PATH(tx)
        }
        
        # Set justifications if changed.
        if (PSI_HJUSTIFY(g_kt) != TX_HJUSTIFY(tx)) {
            call sprintf (g_output, SZ_LINE, "%d HJ")
            call pargi (TX_HJUSTIFY(tx))
            call psk_out (g_output)
            PSI_HJUSTIFY(g_kt) = TX_HJUSTIFY(tx)
        }
        
        if (PSI_VJUSTIFY(g_kt) != TX_VJUSTIFY(tx)) {
            call sprintf (g_output, SZ_LINE, "%d VJ")
            call pargi (TX_VJUSTIFY(tx))
            call psk_out (g_output)
            PSI_VJUSTIFY(g_kt) = TX_VJUSTIFY(tx)
        }
        
        # Set mono-spaced if changed.
        if (PSI_VARIABLE(g_kt) != TX_VARIABLE(tx)) {
            call sprintf (g_output, SZ_LINE, "%b VT")
            call pargb (itob (TX_VARIABLE(tx)))
            call psk_out (g_output)
            PSI_VARIABLE(g_kt) = TX_VARIABLE(tx)
        }
        
        # Set other attributes.
        call psi_color (TX_COLOR(tx))
        call psi_ltype (GL_SOLID)
        
        # Move to the point specified for the text.
        call psi_move (xc, yc)
        
        # Since font changes can occur in the middle of text, process the text
        # in sections, each section corresponding to the current font.
        call salloc (out_text, 3 * n, TY_CHAR)
        call salloc (ps_text, 9 * n, TY_CHAR)
        nchar = stx_segment (text, n, Memc[out_text], TX_FONT(tx))
        
        # The PostScript handles the different paths by reversing the order
        # of characters and just dealing with path RIGHT and UP.
        if (TX_PATH(tx) == GT_LEFT || TX_PATH(tx) == GT_DOWN) {
            call salloc (out_text2, nchar, TY_CHAR)
            call psi_reverse (Memc[out_text], Memc[out_text2])
            out_text = out_text2
        }
        
        # Write the text out.  The text starts with a PostScript mark.
        call psk_out ("MK")
        ip = out_text
        while (Memc[ip] != EOS) {
            
            # Write the font out.
            font = Memc[ip]

            ip = ip + 1
            switch (font) {
            case GT_ROMAN:
                call psk_out ("RF")
            case GT_GREEK:
                call psk_out ("GF")
            case GT_ITALIC:
                call psk_out ("IF")
            case GT_BOLD:
                call psk_out ("BF")
	    case PS_USER:
		call psk_out ("PF")
	    case PS_SETUP:
		call psk_out (Memc[ip])
            default:
                call psk_out ("RF")
            }
            
	# Check to make sure the font is not set to setup a new font.
	# If it is, do not print anything else out, since the setup
	# string has already been printed out above.
	    if (font < PS_SETUP) {
		# Write the text out.
            	call psi_mkpsstring (Memc[ip], Memc[ps_text], 9 * n )
            	call psk_out (Memc[ps_text])
            }

            # Set to the next string.
            while (Memc[ip] != EOS)
                ip = ip + 1
            ip = ip + 1
            
        }
       
        # Now put the PostScript command to draw the string.
        call psk_out ("WS")

        # Set the last font defined.
        PSI_TXFONT(g_kt) = font
        
        # Free memory.
        call sfree (sp)
        
end
#---------------------------------------------------------------------------
# End of psi_text
#---------------------------------------------------------------------------

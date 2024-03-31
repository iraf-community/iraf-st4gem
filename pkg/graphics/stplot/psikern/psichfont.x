include <gset.h>
include "psi.h"

#---------------------------------------------------------------------------
.help psi_change_font 1May92 st4gem.stplot
.ih
NAME
psi_change_font -- Change the PostScript fonts used for the GKI fonts.
.ih
USAGE
call psi_change_font (font, font_name)
.ih
ARGUMENTS
.ls font (int)
Which font to change.  This is one of the GKI codes for each of the
GKI fonts.
.le
.ls font_name (char[ARB])
The name of the PostScript font to use.
.le
.ih
DESCRIPTION
This routine changes the association between GIO font and PostScript font.
A pre-defined association is defined in the graphcap entry for the PSIKern
devices.  Normally, the GIO ROMAN font corresponds to the PostScript font
Times-Roman, etc.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------
procedure psi_change_font (font, font_name)

int     font                    # The GKI font to change.
char    font_name[ARB]          # The PostScript font name.

# Declarations.
include "psi.com"

begin
	# Set DCS page indicator.
	if (g_bpage) {
	    call psk_page
	    g_bpage = false
	}
	
        # Define the new font.
        switch (font) {
        case GT_ROMAN:
            call psk_defl ("RF", font_name)
        case GT_GREEK:
            call psk_defl ("GF", font_name)
        case GT_ITALIC:
            call psk_defl ("IF", font_name)
        case GT_BOLD:
            call psk_defl ("BF", font_name)
        }
        
        # If the association that has changed is for the current font,
        # null the current font.
        if (font == PSI_TXFONT(g_kt))
            PSI_TXFONT(g_kt) = -1
        
        # Add the font to the fontlist.
        call psk_addfont (font_name, PSI_FONTLIST(g_kt))
        
end
#---------------------------------------------------------------------------
# End of psi_change_font
#---------------------------------------------------------------------------

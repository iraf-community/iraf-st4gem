include <gki.h>
include <gset.h>
include "psi.h"

# Define the rotations for Portrait or Landscape.
define PORTRAIT_ROTATION  0.
define LANDSCAPE_ROTATION 90.

# Define conversion from meters to points.
define METER_TO_POINT 72.0 * 39.37

#---------------------------------------------------------------------------
.help psi_prolog 1May92 plot
.ih
NAME
psi_prolog -- Write out the necessary PostScript prolog.
.ih
USAGE
call psi_prolog
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_prolog

# Declarations.
int     lx, ly, ux, uy           # Lower/upper corners of the Bounding Box

include "psi.com"
include "psiparams.com"

# Function Prototypes
long    clktime()
int     strlen()

begin
        
        # Write out the header comments.
        call psk_out ("%!PS-Adobe-1.0")
        call psk_flush
        
        call psk_out ("%%Creator: PSIKern - An IRAF GKI Translator")
        call psk_flush
        
        call strcpy ("%%CreationDate: ", g_output, SZ_LINE)
        call cnvtime (clktime (0), g_output[strlen(g_output)+1], SZ_LINE)
        call psk_out (g_output)
        call psk_flush
        
        call psk_out ("%%DocumentFonts: (atend)")
        call psk_flush
        
        call psk_out ("%%Pages: (atend)")
        call psk_flush
        
        if (PSI_PORTRAIT(g_kt)) {
            lx = METER_TO_POINT * PSI_XOFF(g_kt)
            ly = METER_TO_POINT * PSI_YOFF(g_kt)
            ux = lx +  (METER_TO_POINT * PSI_XSIZE(g_kt))
            uy = ly +  (METER_TO_POINT * PSI_YSIZE(g_kt))
        } else {
            lx = METER_TO_POINT * PSI_YOFF(g_kt)
            ly = METER_TO_POINT * PSI_XOFF(g_kt)
            ux = lx +  (METER_TO_POINT * PSI_YSIZE(g_kt))
            uy = ly +  (METER_TO_POINT * PSI_XSIZE(g_kt))
        }
        call sprintf (g_output, SZ_LINE, "%%%%BoundingBox: %d %d %d %d")
        call pargi (lx)
        call pargi (ly)
        call pargi (ux)
        call pargi (uy)
        call psk_out (g_output)
        call psk_flush
        
        call psk_out ("%%EndComments")
        call psk_flush
        
        # Define page size/orientation.
        call psk_defi ("MAXNDC", PS_OUT_RESOLUTION)
        call psk_defr ("PageSizeXMeter", PSI_XSIZE(g_kt))
        call psk_defr ("PageSizeYMeter", PSI_YSIZE(g_kt))
        call psk_defr ("PageOffXMeter", PSI_XOFF(g_kt))
        call psk_defr ("PageOffYMeter", PSI_YOFF(g_kt))
        call psk_defb ("PortraitMode", PSI_PORTRAIT(g_kt))
        call psk_defr ("PortraitRotation", PORTRAIT_ROTATION)
        call psk_defr ("LandScapeRotation", LANDSCAPE_ROTATION)
        
        # Define some GKI constants to be used in the PostScript code.
        call psk_defi ("GT_LEFT", GT_LEFT)
        call psk_defi ("GT_RIGHT", GT_RIGHT)
        call psk_defi ("GT_TOP", GT_TOP)
        call psk_defi ("GT_BOTTOM", GT_BOTTOM)
        call psk_defi ("GT_UP", GT_UP)
        call psk_defi ("GT_DOWN", GT_DOWN)
        call psk_defi ("GT_CENTER", GT_CENTER)
        
        # Define the fonts.
        PSI_FONTLIST(g_kt) = NULL
        call psi_change_font (GT_ROMAN, roman_font)
        call psi_change_font (GT_GREEK, greek_font)
        call psi_change_font (GT_ITALIC, italic_font)
        call psi_change_font (GT_BOLD, bold_font)
        
        # Define the dot/dash sizes.
        call psk_defi ("DASH", PSI_DASH(g_kt))
        call psk_defi ("DOT", PSI_DOT(g_kt))
        call psk_defi ("SPACE", PSI_SPACE(g_kt))

	# Define the default background color
	if (!IS_INDEFI(PSI_BACK(g_kt)))
	    call psk_defi("BACKCOLOR", PSI_BACK(g_kt))
        
        # End the prolog.
        call psk_flush
        call psk_prolog (PSI_PROLOG(g_kt))
	call psk_flush

	# Set the graphics/image LUT's
	g_maxgrcolor = DEF_MAX_COLOR
	call psi_set_lut (PS_GRAPHIC, lut_gr)
	call psi_set_lut (PS_IMAGE, lut_im)
	g_defmaxgrcolor = g_maxgrcolor
	
	# Set the page scaling and save VM.
	call psk_out ("SS S")

	# End the prolog.
	call psk_flush
	call psk_out ("%%EndProlog")
	call psk_flush
end
#---------------------------------------------------------------------------
# End of psi_prolog
#---------------------------------------------------------------------------

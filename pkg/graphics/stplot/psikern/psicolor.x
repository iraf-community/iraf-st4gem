include "psi.h"

#---------------------------------------------------------------------------
.help psi_color 1Nov94 source
.ih
NAME
psi_color -- Set the graphics drawing color in PostScript.
.ih
USAGE
call psi_color (index)
.ih
ARGUMENTS
.ls index (int)
I: The index into the graphics color table to set the appropriate color.
.le
.ih
DESCRIPTION
The PostScript defines three arrays, GR, GG, and GB, containing,
respectively, the red, green, and blue components for the possible graphics
colors.  The default is the same as that defined by tvmark:
.nf
         0 - black
         1 - white
         2 - red
         3 - green
         4 - blue
         5 - yellow
         6 - cyan
         7 - magenta
         8 - coral
         9 - maroon
        10 - orange
        11 - khaki
        12 - orchid
        13 - turquoise
        14 - violet
        15 - wheat
.fi
The three arrays are read for the red, green, and blue components and a call
to setrgbcolor is made.  This is available in all implementations of 
PostScript, though most printers will translate this to gray scale.
.ih
SEE ALSO
t_psikern, tvmark
.endhelp
#---------------------------------------------------------------------------

procedure psi_color (index)

int index  # I:  The color index.

# Declarations
include "psi.com"

begin
	if (g_debug)
	    call eprintf ("psi_color: Setting color.\n")
	
	# Set DCS page indicator.
	if (g_bpage) {
	    call psk_page
	    g_bpage = false
	}
	
        # Check for valid color.
        if (index > g_maxgrcolor) {
	    if (g_debug) {
		call eprintf ("psikern: Color %d greater than maximum %d\n")
		call pargi (index)
		call pargi (g_maxgrcolor)
	    }
            index = DEF_COLOR
        }
        if (index < 0) {
	    if (g_debug) {
		call eprintf ("psikern: Color %d is less than 0, setting to default\n")
		call pargi (index)
	    }
            index = DEF_COLOR
        }

        # Check to see if the color has changed.  If so, tell the device.
        if (PSI_COLOR (g_kt) != index) {

            call sprintf (g_output, SZ_LINE, "%d SC")
            call pargi (index)
            call psk_out (g_output)
            
            PSI_COLOR(g_kt) = index
        }
end
#---------------------------------------------------------------------------
# End of psi_color
#---------------------------------------------------------------------------

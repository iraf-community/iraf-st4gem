include	<error.h>
include	<psiescape.h>
include	"psi.h"

# Memory management.
define	R			Mems[r]
define	G			Mems[g]
define 	B			Mems[b]

#---------------------------------------------------------------------------
.help psi_set_lut 4Nov94 source
.ih
NAME
psi_set_lut - Set LUT from a table.
.endhelp
#---------------------------------------------------------------------------
procedure psi_set_lut (type, table)

int	type			# I:  Type of LUT to set.
char	table[ARB]		# I:  Table to use.

# Declarations
include	"psi.com"

int	n			# Number of colors.
pointer	r, g, b, name		# red/green/blue/name of colors.
int	strlen()		# Length of string.

begin
	# Make sure a table is specified.  If not, nothing will be done.
	if (strlen(table) > 0) {
	    
	    # Read in the colormap from the table.
	    iferr (call psi_lut_read (table, r, g, b, name, n)) {
		call eprintf ("psikern: could not read colormap table %s\n")
		call pargstr (table)
		call erract (EA_WARN)
	    }

	    # Else, set the colormaps.
	    else {
		
		# Set the graphics colormap.
		if (type == PS_GRAPHIC) {
		    call psi_download_lut (PS_GR_RED_LUT, R, n)
		    call psi_download_lut (PS_GR_GREEN_LUT, G, n)
		    call psi_download_lut (PS_GR_BLUE_LUT, B, n)
		}

		# Set the image colormap.
		else {
		    call psi_download_lut (PS_IMAGE_RED_LUT, R, n)
		    call psi_download_lut (PS_IMAGE_GREEN_LUT, G, n)
		    call psi_download_lut (PS_IMAGE_BLUE_LUT, B, n)
		}

		# Free the arrays.
		call mfree (r, TY_SHORT)
		call mfree (g, TY_SHORT)
		call mfree (b, TY_SHORT)
		call mfree (name, TY_CHAR)
	    }
	}
end
#---------------------------------------------------------------------------
# End of psi_set_lut
#---------------------------------------------------------------------------

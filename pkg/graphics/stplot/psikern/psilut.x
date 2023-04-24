include <psiescape.h>
include "psi.h"

define MAX_CHAR 72

#---------------------------------------------------------------------------
.help psi_download_lut 1May92 plot
.ih
NAME
psi_download_lut -- Set lookup tables for images and graphics.
.ih
USAGE
call psi_download_lut (lut_type, lut, lut_length)
.ih
ARGUMENTS
.ls lut_type (int)
Which lookup table to set.  This is equavelent to the function code in the
GKI escape sequence.  The definitions can be found in the include file
psiescape.h
.le
.ls lut (short[length_lut])
The new lookup table.  The contents of the lookup table is an integer from
the range 0 (no color) to 255 (full color).  Since LUTs are usually specified
from 0-1, there are two macros defined in psiescape.h to help in scaling from
the normalized system to the GKI system call PS_PACKLUT and PS_UNPACKLUT.
See psiescape.h for more info.
.le
.ls lut_length (int)
The number of values in the lookup table.
.le
.ih
DESCRIPTION
This routine is called to define/redefine the image and graphics lookup tables
in the PostScript code to change the way graphics and images appear in the
output.  There is a default graphics LUT which is that defined by the IRAF task
tvmark.  For images, no default is defined.  If an image is displayed without
defining a lookup table, a grey-scale transformation defined within PostScript
is used, else the user-defined lookup table is used.

The size of the graphics LUT's really don't matter much, as long as all three
color components are the same length and at least contain PS_GR_LUT_SIZE
entries.  However, the image LUT's should be PS_IMAGE_LUT_SIZE exactly.
The PostScript color image operators wouldn't be able to deal with anything
else.
.ih
SEE ALSO
t_psikern
.endhelp
#---------------------------------------------------------------------------

procedure psi_download_lut (lut_type, lut, lut_length)

int     lut_type                # Type of LUT to set.
short   lut[ARB]	        # The LUT.
int     lut_length              # The length of the LUT.

# Declarations.
int     i, j                    # More counters.
int     index                   # Index into the LUT.
int     ip                      # Pointer into output string.
int     nlines                  # Number of lines of table values.
int     nvalues                 # Number of table values per line.

include "psi.com"

begin
        # If this is for the graphics color palette, retrieve the
        # size of the lookup table.
        switch (lut_type){
        case PS_GR_RED_LUT, PS_GR_GREEN_LUT, PS_GR_BLUE_LUT:
            g_maxgrcolor = max (0, lut_length - 1)
	    PSI_COLOR(g_kt) = -1
        }
        
        # Write out the appropraite name.
        switch (lut_type) {
        case PS_IMAGE_RED_LUT:
            call psk_out ("/IR")
        case PS_IMAGE_GREEN_LUT:
            call psk_out ("/IG")
        case PS_IMAGE_BLUE_LUT:
            call psk_out ("/IB")
        case PS_GR_RED_LUT:
            call psk_out ("/GR")
        case PS_GR_GREEN_LUT:
            call psk_out ("/GG")
        case PS_GR_BLUE_LUT:
            call psk_out ("/GB")
        }
        call sprintf (g_output, SZ_LINE, "%d MI")
        call pargi (lut_length)
        call psk_out (g_output)
        
        # Give it the data.
        nvalues = MAX_CHAR / 2
        nlines = lut_length / nvalues
        index = 1
        ip = 1
        do i = 1, nlines {
            do j = 0, nvalues - 1 {
                call sprintf (g_output[ip], MAX_CHAR, "%02.2x")
                call pargs (lut[index + j])
                ip = ip + 2
            }
            call psk_out (g_output)
            ip = 1
            index = index + nvalues
        }
        ip = 1
        for (i = index; i <= lut_length; i = i + 1) {
            call sprintf (g_output[ip], MAX_CHAR, "%02.2x")
            call pargs (lut[i])
            ip = ip + 2
        }
        call psk_out (g_output)
        
end
#---------------------------------------------------------------------------
# End of psi_download_lut
#---------------------------------------------------------------------------

include "psi.h"

# Define the number of bits per sample the PostScript image operator
# will require.
define BITSPERSAMPLE 8

define MAX_CHAR 72

#---------------------------------------------------------------------------
.help psi_putcellarray 1May92 plot
.NAME
psi_putcellarray -- Draw the cell array out to the PostScript device.
.ih
USAGE
call psi_putcellarray (m, nx, ny, ax1, ay1, ax2, ay2) 
.ih
ARGUMENTS
.ls m (short[nx,ny])
I: The cell array.
.le
.ls nx, ny (int)
I: The number of pixels in the X, Y dimension
.le
.ls ax1,ay1 (int)
I: The lower left corner of the output window.
.le
.ls ax2, ay2 (int)
I: The upper, right corner of the output window.
.le
.ih
DESCRIPTION
This routine mapes the 8bit image into PostScript that is then imaged onto
the page using the IMAGE operator.
.ih
SEE ALSO
t_psikern
.endhelp
#
# 10-18-96 Perry Greenfield	Fixed problem caused by previous fix. Caused
#				small images to fail.
#  9-11-96 Jon Eisenhamer	Fixed problem when image lines were multiple
#				of 36 pixels.
#  3-13-97 Warren Hack		Revised section that converted short's into
#				ASCII hex by replacing 'sprintf' calls on
#				each pixel by string lookup replacement.
#				This resulted in a 10x speed increase.
#---------------------------------------------------------------------------

procedure psi_putcellarray (m, nx, ny, ax1, ay1, ax2, ay2)

short m[ARB]                    # I:  The cell array.
int   nx, ny                    # I:  Size in X, Y of the cell array.
int   ax1, ay1                  # I:  The lower, left corner of the output 
                                #     window.
int   ax2, ay2                  # I:  The upper, right corner of the output 
                                #     window.

# Declarations
real dx, dy                     # Change in the window viewport.
real sax1, sax2, say1, say2     # Rescaled input coordinates.
real sx, sy                     # The scale factor between the window and 
                                # image.
real tx, ty                     # The translation factor between the window 
                                # and image.

int i, j                        # Indexes into the image cell array.
int index                       # The index into the image cell array.
int ip                          # Index inout output buffer.
int nlines                      # Number of full lines of image cell values.
int nvalues                     # Number of image cell values that can be on 
                                # a line.
int total                       # Total number of pixels in the image.

string hexval "000102030405060708090A0B0C0D0E0F101112131415161718191A1B1C1D1E1F\
202122232425262728292A2B2C2D2E2F303132333435363738393A3B3C3D3E3F\
404142434445464748494A4B4C4D4E4F505152535455565758595A5B5C5D5E5F\
606162636465666768696A6B6C6D6E6F707172737475767778797A7B7C7D7E7F\
808182838485868788898A8B8C8D8E8F909192939495969798999A9B9C9D9E9F\
A0A1A2A3A4A5A6A7A8A9AAABACADAEAFB0B1B2B3B4B5B6B7B8B9BABBBCBDBEBF\
C0C1C2C3C4C5C6C7C8C9CACBCCCDCECFD0D1D2D3D4D5D6D7D8D9DADBDCDDDEDF\
E0E1E2E3E4E5E6E7E8E9EAEBECEDEEEFF0F1F2F3F4F5F6F7F8F9FAFBFCFDFEFF"

string endchar " R"

int	chridx
include "psi.com"

# Function prototypes
bool fp_equalr()

begin
        # If this feature is not supported, don't do it.
	if (PSI_DOCELL(gkt) == NO) {
	    if (g_debug)
		call eprintf ("psi_putcellarray: Not done because 'ca' is not set\n")
	    return
	}
	
        if (g_debug){
            call eprintf ("psi_putcellarray: Array size %d, %d, at %d, %d to %d, %d.\n")
            call pargi (nx)
            call pargi (ny)
            call pargi (ax1)
            call pargi (ay1)
            call pargi (ax2)
            call pargi (ay2)
        }
        
        # Count as a drawing operation.
        g_ndraw = g_ndraw + 1

	# Set DCS page indicator.
	if (g_bpage) {
	    call psk_page
	    g_bpage = false
	}
	
        # Rescale the input coordinates.
        sax1 = ax1 * PSI_GKI2OUT(g_kt)
        sax2 = ax2 * PSI_GKI2OUT(g_kt)
        say1 = ay1 * PSI_GKI2OUT(g_kt)
        say2 = ay2 * PSI_GKI2OUT(g_kt)
        
        # Determine the image transformation matrix.
        dx = sax2 - sax1
        dy = say2 - say1
        if (fp_equalr (dx, 0.0))
            sx = 1.
        else
            sx = nx / dx
        
        if (fp_equalr (dy, 0.0))
            sy = 1.
        else
            sy = ny / dy
        
        tx = -(sx * sax1)
        ty = -(sy * say1)
        
        # Start the cell drawing command.
        call sprintf (g_output, MAX_CHAR, "S %d %d %d [%g 0 0 %g %g %g] PC")
        call pargi (nx)
        call pargi (ny)
        call pargi (BITSPERSAMPLE)
        call pargr (sx)
        call pargr (sy)
        call pargr (tx)
        call pargr (ty)
        call psk_out (g_output)
        
        # Give it the data.
        nvalues = MAX_CHAR / 2
        total = nx * ny
        nlines = (total-1) / nvalues
        index = 1
        ip = 1
        do i = 1, nlines {
            do j = 0, nvalues - 1 {
                #call sprintf (g_output[ip], MAX_CHAR, "%02.2x")
                #call pargs (m[index + j])
		chridx = (m[index + j] * 2) + 1 
		g_output[ip] = hexval[chridx]
		g_output[ip+1] = hexval[chridx+1]
                ip = ip + 2
            }
            call psk_out (g_output)
            ip = 1
            index = index + nvalues
        }

        ip = 1
        for (i = index; i <= total; i = i + 1) {
            #call sprintf (g_output[ip], MAX_CHAR, "%02.2x")
            #call pargs (m[i])
	    chridx = (m[i] * 2) + 1
	    g_output[ip] = hexval[chridx]
	    g_output[ip+1] = hexval[chridx+1]
            ip = ip + 2
        }
	#call sprintf(g_output[ip],MAX_CHAR," R")
	g_output[ip] = endchar[1]
	g_output[ip+1] = endchar[2]
	g_output[ip+2] = EOS
        call psk_out (g_output)
        #call psk_out ("R")
        
end
#---------------------------------------------------------------------------
# End of psi_putcellarray
#---------------------------------------------------------------------------

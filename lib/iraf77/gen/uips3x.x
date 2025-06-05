include <imhdr.h>
include <iraf77.h>

# UIPS3? -- put a section to an apparently 3-d file

procedure uips3x (im, xmin, xmax, ymin, ymax, zmin, zmax, buffer, istat)

pointer im           # pointer to the image header file
int xmin             # x-pixel lower limit
int xmax             # x-pixel upper limit
int ymin             # y-pixel lower limit
int ymax             # y-pixel upper limit
int zmin             # z-pixel lower limit
int zmax             # z-pixel upper limit
complex buffer[ARB]    # user's buffer to be transfered to imio buffer
int istat            # return status code

pointer ip, imps3x()
int npix

begin
	istat = ER_OK

	# get the size of a section
	npix = xmax - xmin + 1
	npix = npix * (ymax - ymin + 1)
	npix = npix * (zmax - zmin + 1)

	# call imio to get an imio buffer area
	iferr (ip = imps3x(im, xmin, xmax, ymin, ymax, zmin, zmax))
		{istat = ER_IMWRITE
		return
		}

	# load the user's array into imio buffer
	call amovx (buffer, Memx[ip], npix)

end

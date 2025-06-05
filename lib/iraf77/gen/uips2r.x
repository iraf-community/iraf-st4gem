include <imhdr.h>
include <iraf77.h>

# UIPS2? -- put a section to an apparently 2-d file

procedure uips2r (im, xmin, xmax, ymin, ymax, buffer, istat)

pointer im           # pointer to the image header file
int xmin             # x-pixel lower limit
int xmax             # x-pixel upper limit
int ymin             # y-pixel lower limit
int ymax             # y-pixel upper limit
real buffer[ARB]    # user's buffer to be transfered to imio buffer
int istat            # return status code

pointer ip, imps2r()
int npix

begin
	istat = ER_OK

	# get the size of a section
	npix = xmax - xmin + 1
	npix = npix * (ymax - ymin + 1)

	# call imio to get an imio buffer area
	iferr (ip = imps2r(im, xmin, xmax, ymin, ymax))
		{istat = ER_IMWRITE
		return
		}

	# load the user's array into imio buffer
	call amovr (buffer, Memr[ip], npix)

end

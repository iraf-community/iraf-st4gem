include <imhdr.h>
include <iraf77.h>

# UIGS2? -- get a section from an apparently 2-d file

procedure uigs2l (im, xmin, xmax, ymin, ymax, buffer, istat)

pointer im           # pointer to the image header file
int xmin             # x-pixel lower limit
int xmax             # x-pixel upper limit
int ymin             # y-pixel lower limit
int ymax             # y-pixel upper limit
long buffer[ARB]    # user's buffer to be filled by this routine
int istat            # return status code

pointer ip, imgs2l()
int npix

begin
	istat = ER_OK

	# get the size of a section
	npix = xmax - xmin + 1
	npix = npix * (ymax - ymin + 1)

	# call imio to read the data into an imio buffer area
	iferr (ip = imgs2l(im, xmin, xmax, ymin, ymax)) {
		istat = ER_IMREAD
		return
		}
	# load the user's array from imio buffer
	call amovl (Meml[ip], buffer, npix)

end

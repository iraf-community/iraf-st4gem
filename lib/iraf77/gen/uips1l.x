include <imhdr.h>
include <iraf77.h>

# UIPS1? -- put a data line subsection to an apparently 1-d file

procedure uips1l (im, xmin, xmax, buffer, istat)

pointer im           # pointer to the image header file
int xmin             # x-pixel lower limit
int xmax             # x-pixel upper limit
long buffer[ARB]    # user's buffer to be transfered to imio buffer
int istat            # return status code

pointer ip, imps1l()
int npix

begin
	istat = ER_OK

	# get the size of a line subsection
	npix = xmax - xmin + 1

	# call imio to get an imio buffer area
	iferr (ip = imps1l(im, xmin, xmax))
		{istat = ER_IMWRITE
		return
		}

	# load the user's array into imio buffer
	call amovl (buffer, Meml[ip], npix)

end

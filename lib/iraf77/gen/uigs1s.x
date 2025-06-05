include <imhdr.h>
include <iraf77.h>

# UIGS1? -- get a data line subsection from an apparently 1-d file

procedure uigs1s (im, xmin, xmax, buffer, istat)

pointer im           # pointer to the image header file
int xmin             # x-pixel lower limit
int xmax             # x-pixel upper limit
short buffer[ARB]    # user's buffer to be filled by this routine
int istat            # return status code

pointer ip, imgs1s()
int npix

begin
	istat = ER_OK

	# get the size of a line subsection
	npix = xmax - xmin + 1

	# call imio to read the data into an imio buffer area
	iferr (ip = imgs1s(im, xmin, xmax))
		{istat = ER_IMREAD
		return
		}

	# load the user's array from imio buffer
	call amovs (Mems[ip], buffer, npix)

end

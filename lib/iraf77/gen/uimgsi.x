include <imhdr.h>
include <iraf77.h>

# UIMGS? -- get a section from a n-dimensional file(greater than 2)

procedure uimgsi (im, fpix, lpix, buffer, istat)

pointer im           # pointer to the image header file
int fpix[7]          # array giving starting pixels to get from image
int lpix[7]          # array giving last pixels to get from image
int buffer[ARB]    # user's buffer to be filled by this routine
int istat            # return status code

pointer ip, imggsi()
int npix, ndim
long vs[IM_MAXDIM], ve[IM_MAXDIM]
int i

begin
	istat = ER_OK

	# get the dimensionality of the image
	ndim = IM_NDIM(im)

	# set up the starting and ending vectors
	# and get the size of a section
	npix = 1
	do i = 1, ndim {
		vs[i] = fpix[i]
		ve[i] = lpix[i]
                npix = npix * (lpix[i] - fpix[i] + 1)
	}

	# call imio to read the data into an imio buffer area
	iferr (ip = imggsi(im, vs, ve, ndim)) {
		istat = ER_IMREAD
		return
		}
	# load the user's array from imio buffer
	call amovi (Memi[ip], buffer, npix)

end

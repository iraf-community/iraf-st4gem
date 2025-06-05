include <imhdr.h>
include <iraf77.h>

# UIMPS? -- put a section to a n-dimensional file(greater than 2)

procedure uimpsi (im, fpix, lpix, buffer, istat)

pointer im           # pointer to the image header file
int fpix[7]          # array giving starting pixels to get from image
int lpix[7]          # array giving last pixels to get from image
int buffer[ARB]    # user's buffer to be transfered to imio buffer
int istat            # return status code

pointer ip, impgsi()
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

	# call imio to get an imio buffer area
	iferr (ip = impgsi(im, vs, ve, ndim))
		{istat = ER_IMWRITE
		return
		}
	# load the user's array into imio buffer
	call amovi (buffer, Memi[ip], npix)

end

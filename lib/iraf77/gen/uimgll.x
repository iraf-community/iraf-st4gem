include <imhdr.h>
include <iraf77.h>

# UIMGL? -- get a data line from a n-dimensional file(greater than 2)

procedure uimgll (im, linenum, buffer, istat)

pointer im           # pointer to the image header file
int linenum[6]       # array describing line number to get from image
long buffer[ARB]    # user's buffer to be filled by this routine
int istat            # return status code

pointer ip, imggsl()
int ndim
long vs[IM_MAXDIM], ve[IM_MAXDIM]
int i

begin
	istat = ER_OK

	# get the dimensionality of the image
	ndim = IM_NDIM(im)

	# set up the starting and ending vectors
	vs[1] = 1
	ve[1] = IM_LEN(im, 1)

	do i = 2, ndim {
		vs[i] = linenum[i-1]
		ve[i] = linenum[i-1]
	}

	# call imio to read the data into an imio buffer area
	iferr (ip = imggsl(im, vs, ve, ndim)) {
		istat = ER_IMREAD
		return
		}
	# load the user's array from imio buffer
	call amovl (Meml[ip], buffer, IM_LEN(im,1))

end

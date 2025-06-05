include <imhdr.h>
include <iraf77.h>

# UIMPL? -- put a data line to a n-dimensional file(greater than 2)

procedure uimpll (im, linenum, buffer, istat)

pointer im           # pointer to the image header file
int linenum[6]       # array describing line number to put to image
long buffer[ARB]    # user's buffer to be transfered to imio buffer
int istat            # return status code

pointer ip, impgsl()
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

	# call imio to get an imio buffer area
	iferr (ip = impgsl(im, vs, ve, ndim))
		{istat = ER_IMWRITE
		return
		}
	# load the user's array into imio buffer
	call amovl (buffer, Meml[ip], IM_LEN(im,1))

end

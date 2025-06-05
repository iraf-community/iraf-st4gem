include <imhdr.h>
include <iraf77.h>

# UIPL1? -- put a data line to an apparently 1-d file

procedure uipl1r (im, buffer, istat)

pointer im           # pointer to the image header file
real buffer[ARB]    # user's buffer to be transfered to imio buffer
int istat            # return status code

pointer ip, impl1r()

begin
	istat = ER_OK

	# call imio to get an imio buffer area
	iferr (ip = impl1r(im))
		{istat = ER_IMWRITE
		return
		}

	# load the user's array into imio buffer
	call amovr (buffer, Memr[ip], IM_LEN(im,1))

end

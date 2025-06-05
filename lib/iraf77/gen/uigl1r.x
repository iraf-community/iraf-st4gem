include <imhdr.h>
include <iraf77.h>

# UIGL1? -- get a data line from an apparently 1-d file

procedure uigl1r (im, buffer, istat)

pointer im           # pointer to the image header file
real buffer[ARB]    # user's buffer to be filled by this routine
int istat            # return status code

pointer ip, imgl1r()

begin
	istat = ER_OK

	# call imio to read the data into an imio buffer area
	iferr (ip = imgl1r(im))
		{istat = ER_IMREAD
		return
		}

	# load the user's array from imio buffer
	call amovr (Memr[ip], buffer, IM_LEN(im,1))

end

include <imhdr.h>
include <iraf77.h>

# UIGL2? -- get a data line from an apparently 2-d file

procedure uigl2r (im, linenum, buffer, istat)

pointer im           # pointer to the image header file
int linenum          # line number to get from image
real buffer[ARB]    # user's buffer to be filled by this routine
int istat            # return status code

pointer ip, imgl2r()

begin
	istat = ER_OK

	# call imio to read the data into an imio buffer area
	iferr (ip = imgl2r(im, linenum)) {
		istat = ER_IMREAD
		return
		}
	# load the user's array from imio buffer
	call amovr (Memr[ip], buffer, IM_LEN(im,1))

end

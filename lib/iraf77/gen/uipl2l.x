include <imhdr.h>
include <iraf77.h>

# UIPL2? -- put a data line to an apparently 2-d file

procedure uipl2l (im, linenum, buffer, istat)

pointer im           # pointer to the image header file
int linenum          # line number to put to image
long buffer[ARB]    # user's buffer to be transfered to imio buffer
int istat            # return status code

pointer ip, impl2l()

begin
	istat = ER_OK

	# call imio to get an imio buffer area
	iferr (ip = impl2l(im, linenum))
		{istat = ER_IMWRITE
		return
		}
	# load the user's array into imio buffer
	call amovl (buffer, Meml[ip], IM_LEN(IM,1))

end

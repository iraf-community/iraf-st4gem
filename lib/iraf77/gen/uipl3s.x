include <imhdr.h>
include <iraf77.h>

# UIPL3? -- put a data line to an apparently 3-d file

procedure uipl3s (im, linenum, bandnum, buffer, istat)

pointer im           # pointer to the image header file
int linenum          # line number to put to image
int bandnum          # band number to put to image
short buffer[ARB]    # user's buffer to be transfered to imio buffer
int istat            # return status code

pointer ip, impl3s()

begin
	istat = ER_OK

	# call imio to get an imio buffer area
	iferr (ip = impl3s(im, linenum, bandnum))
		{istat = ER_IMWRITE
		return
		}

	# load the user's array into imio buffer
	call amovs (buffer, Mems[ip], IM_LEN(im,1))

end

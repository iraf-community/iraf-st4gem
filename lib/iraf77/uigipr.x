include <imio.h>
include <imset.h>
include <iraf77.h>

# UIGIPR--Get image processing parameter of type real

procedure uigipr (im, parnum, buffer, istat)

pointer im	# image descriptor pointer returned by uimopn
int parnum	# symbolic code; image parameter to get
real buffer	# buffer for parameter value
int istat	# return status code

# For the present, only the out-of-bounds pixel constant in imio can be
# set or gotten by the user in imio.

begin

	switch (parnum){
	case IM_BNDRYPIXVAL:
		buffer = IM_OOBPIX(im)
		istat = ER_OK

	default:
		buffer = 0.
		istat = ER_IMSTATUNKPAR
	}

	return

end

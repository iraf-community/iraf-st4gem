include <imio.h>
include <imset.h>
include <iraf77.h>

# UIPIPR--Set image processing parameter of type real

procedure uipipr (im, parnum, buffer, istat)

pointer im	# image descriptor pointer returned by uimopn
int parnum	# symbolic code; image parameter to get
real buffer	# buffer for parameter value
int istat	# return status code

# For the present, only the out-of-bounds pixel constant in imio can be
# set or gotten by the user in imio.

begin

	istat = ER_OK

	switch (parnum){
	case IM_BNDRYPIXVAL:
	    iferr (call imsetr (im, parnum, buffer))
		istat = ER_IMSETUNKPAR

	default:
	    istat = ER_IMSETUNKPAR
	}

	return

end

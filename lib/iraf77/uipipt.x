include	<imset.h>
include	<imio.h>
include	<iraf77.h>

# UIPIPT -- PUT an IMIO option of type string.

procedure uipipt (im, parnum, buffer, istat)

pointer	im			# image descriptor returned by uimopn
int	parnum			# symbolic code; image parameter to set
%	character*(*) buffer
int	istat			# return status code

int	maxch

# For the present, only the image name in imio can be
# set or gotten by the user in imio.

begin

	switch (parnum) {
	case IM_IMAGENAME:
	    maxch = SZ_IMNAME
	    call f77upk (buffer, IM_NAME(im), maxch)
	    istat = ER_OK
	default:
	    istat = ER_IMSETUNKPAR
	}
end

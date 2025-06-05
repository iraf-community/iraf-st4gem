include <imio.h>
include <imset.h>
include <iraf77.h>

# UIPIPI--Set an image processing parameter of type int

procedure uipipi (im, parnum, buffer, istat)

pointer im	# image descriptor pointer returned by uimopn
int parnum	# symbolic code; image parameter to set
int buffer	# buffer for parameter value
int istat	# return status code

# For the present, only a subset of those parameters declared by imio to be
# set or gotten by the user are allowed by this routine.  The 
# values in the case statement should be compared to the list of parameters
# given in lib$imset.h.  As a specific example, the documentation for 
# this interface originally mentioned fast i/o as one example; imio
# sets this parameter based on INTERNAL criteria so we will not allow it.

begin

	istat = ER_OK

	switch (parnum) {
	case IM_NBNDRYPIX, IM_TYBNDRY, IM_FLAGBADPIX, IM_WHEADER:
	    iferr (call imseti (im, parnum, buffer)) # shouldn't happen!
		istat = ER_IMSETUNKPAR
	default:
	    istat = ER_IMSETUNKPAR
	}

end

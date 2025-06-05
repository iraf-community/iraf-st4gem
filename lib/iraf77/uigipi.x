include <imio.h>
include <imset.h>
include <iraf77.h>

# UIGIPI--Get image processing parameter of type int

procedure uigipi (im, parnum, buffer, istat)

pointer im	# image descriptor pointer returned by uimopn
int parnum	# symbolic code; image parameter to get
int buffer	# buffer for parameter value
int istat	# return status code

int imstati()

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
	    iferr (buffer = imstati (im, parnum)) # shouldn't happen!
		istat = ER_IMSTATUNKPAR
	default:
	    istat = ER_IMSTATUNKPAR
	}

end

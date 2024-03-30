include	"../limit.h"
include "dtrans.h"

#* HISTORY *
#* B.Simon	22-Jul-94	original

# LIMDTRANS -- Compute plot limits for transformation data

procedure limdtrans (dtrans, inlimit, outlimit)

pointer	dtrans		# i: transformation data descriptor
real	inlimit[4]	# i: input plot limits
real	outlimit[4]	# o: calculated plot limits
#--
int	neff, ieff
pointer	xeff, yeff
real	xmin, xmax, ymin, ymax

begin
	if (dtrans == NULL)
	    return

	neff = TRN_NEFF(dtrans)
	xeff = TRN_XEFF(dtrans)
	yeff = TRN_YEFF(dtrans)

	# Get extreme values of the photometric data

	if (neff == 0) {
	    xmin = INDEFR
	    xmax = INDEFR
	    ymin = INDEFR
	    ymax = INDEFR

	} else {
	    xmin = Memr[xeff]
	    xmax = Memr[xeff]
	    ymin = Memr[yeff]
	    ymax = Memr[yeff]

	    do ieff = 2, neff {
		xmin = min (xmin, Memr[xeff+ieff-1])
		xmax = max (xmax, Memr[xeff+ieff-1])
		ymin = min (ymin, Memr[yeff+ieff-1])
		ymax = max (ymax, Memr[yeff+ieff-1])
	    }
	}

	# Update the output limits with these values 
	# where the input values are INDEF

	if (IS_INDEFR(inlimit[LEFT])) {
	    outlimit[LEFT] = xmin
	} else {
	    outlimit[LEFT] = inlimit[LEFT]
	}

	if (IS_INDEFR(inlimit[RIGHT])) {
	    outlimit[RIGHT] = xmax
	} else {
	    outlimit[RIGHT] = inlimit[RIGHT]
	}

	if (IS_INDEFR(inlimit[BOTTOM])) {
	    outlimit[BOTTOM] = ymin
	} else {
	    outlimit[BOTTOM] = inlimit[BOTTOM]
	}

	if (IS_INDEFR(inlimit[TOP])) {
	    outlimit[TOP] = ymax
	} else {
	    outlimit[TOP] = inlimit[TOP]
	}
end

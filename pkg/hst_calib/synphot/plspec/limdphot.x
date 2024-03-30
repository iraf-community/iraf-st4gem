include "dphot.h"
include "../limit.h"

#* HISTORY *
#* B.Simon	08-Jun-94	original
#* B.Simon	14-Jul-94	updated to handle indefs

# LIMDPHOT -- Compute plot limits from photometry data

procedure limdphot (dphot, form, inlimit, outlimit)

pointer	dphot		# u: spectrum descriptor
char	form[ARB]	# i: output form
real	inlimit[4]	# i: input plot limits
real	outlimit[4]	# u: calculated plot limits
#--
int	nwave, istat
pointer	wave, stim
real	minval, maxval

begin
	if (dphot == NULL)
	    return

	# Compute extreme values from photometric statistics

	minval = INDEFR
	maxval = INDEFR
	stim = PHT_STIM(dphot)

	do istat = 0, PHT_NSTAT(dphot) - 1 {
	    if (IS_INDEFR (Memr[stim+istat]))
		next

	    if (IS_INDEFR(minval)) {
		minval = Memr[stim+istat]
	    } else {
		minval = min (minval, Memr[stim+istat])
	    }

	    if (IS_INDEFR(maxval)) {
		maxval = Memr[stim+istat]
	    } else {
		maxval = max (maxval, Memr[stim+istat])
	    }
	}

	if (! IS_INDEFR(minval)) {
	    if (minval == maxval) {
		if (minval > 0.0) {
		    minval = 0.0
		} else {
		    maxval = 0.0
		}
	    }
	}

	# Update horizontal limits from wavelength array

	nwave = PHT_NWAVE(dphot)
	wave = PHT_WAVE(dphot)

	if (IS_INDEFR(inlimit[LEFT])) {
	    if (IS_INDEFR(outlimit[LEFT])) {
		outlimit[LEFT] = Memr[wave]
	    } else {
		outlimit[LEFT] = min (Memr[wave], outlimit[LEFT])
	    }
	} else {
	    outlimit[LEFT] = inlimit[LEFT]
	}

	if (IS_INDEFR(inlimit[RIGHT])) {
	    if (IS_INDEFR(outlimit[RIGHT])) {
		outlimit[RIGHT] = Memr[wave+nwave-1]
	    } else {
		outlimit[RIGHT] = max (Memr[wave+nwave-1], outlimit[RIGHT])
	    }
	} else {
	    outlimit[RIGHT] = inlimit[RIGHT]
	}

	# Update vertical limits from extreme values

	if (IS_INDEFR (inlimit[BOTTOM]) && ! IS_INDEFR(minval)) {
	    if (IS_INDEFR (outlimit[BOTTOM])) {
		outlimit[BOTTOM] = minval
	    } else {
		outlimit[BOTTOM] = min (minval, outlimit[BOTTOM])
	    }
	} else {
	    outlimit[BOTTOM] = inlimit[BOTTOM]
	}

	if (IS_INDEFR (inlimit[TOP]) && ! IS_INDEFR(maxval)) {
	    if (IS_INDEFR (outlimit[TOP])) {
		outlimit[TOP] = maxval
	    } else {
		outlimit[TOP] = max (maxval, outlimit[TOP])
	    }
	} else {
	    outlimit[TOP] = inlimit[TOP]
	}

	# Map the passabd into the range of the plot

	call mapdband (dphot, form, 0.2, outlimit)

end

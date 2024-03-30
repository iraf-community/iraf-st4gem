include "../limit.h"

# SETLIMIT -- Compute the extrema of a passband or spectrum

procedure setlimit (nwave, wave, flux, inlimit, outlimit)

int	nwave		# i: Length of wave and spectrum arrays
real	wave[ARB]	# i: Wavelength array
real	flux[ARB]	# i: Flux array
real	inlimit[4]	# i: Input plot limits
real	outlimit[4]	# u: Output plot limis
#--
int	iwave
real	minval, maxval

begin
	# Update limits from wavelength array

	if (IS_INDEFR(inlimit[LEFT])) {
	    if (IS_INDEFR(outlimit[LEFT])) {
		outlimit[LEFT] = wave[1]
	    } else {
		outlimit[LEFT] = min (wave[1], outlimit[LEFT])
	    }
	} else {
	    outlimit[LEFT] = inlimit[LEFT]
	}

	if (IS_INDEFR(inlimit[RIGHT])) {
	    if (IS_INDEFR(outlimit[RIGHT])) {
		outlimit[RIGHT] = wave[nwave]
	    } else {
		outlimit[RIGHT] = max (wave[nwave], outlimit[RIGHT])
	    }
	} else {
	    outlimit[RIGHT] = inlimit[RIGHT]
	}

	# Compute minimum and maximum of flux array

	minval = INDEFR
	maxval = INDEFR
	do iwave = 1, nwave {
	    if (wave[iwave] < outlimit[LEFT] || wave[iwave] > outlimit[RIGHT])
		next

	    if (IS_INDEFR(flux[iwave]))
		next

	    if (IS_INDEFR(minval)) {
		minval = flux[iwave]
	    } else {
		minval = min (minval, flux[iwave])
	    }

	    if (IS_INDEFR(maxval)) {
		maxval = flux[iwave]
	    } else {
		maxval = max (maxval, flux[iwave])
	    }
	}

	# Update limits from flux array

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

end

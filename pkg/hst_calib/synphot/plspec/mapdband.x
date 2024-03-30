include "dphot.h"
include "../limit.h"

# MAPDBAND -- Map the bandpass function into a new range for plotting
 
procedure mapdband (dphot, form, factor, limit)

pointer	dphot		# u: bandpass descriptor
char	form[ARB]	# i: spectral form
real	factor		# i: proportion of range bandpass will occupy
real	limit[4]	# i: range limits
#--
int	nwave, iwave, iband
real	minval, maxval, maxband, frac
pointer	band

int	is_magunit()
real	ahivr()

begin

	if (dphot == NULL || IS_INDEFR(limit[BOTTOM]) || IS_INDEFR(limit[TOP]))
	    return

	# Compute limits for bandpass

	if (is_magunit (form) == NO) {
	    minval = limit[BOTTOM]
	    maxval = factor * (limit[TOP] - limit[BOTTOM]) + limit[BOTTOM]
	} else {
	    minval = limit[TOP]
	    maxval = (1.0 - factor) * (limit[TOP] - limit[BOTTOM]) + 
		     limit[BOTTOM]
	}

	# Normalize plot between limits

	band = PHT_BAND(dphot)
	nwave = PHT_NWAVE(dphot)
	do iband = 1, PHT_NBAND(dphot) {
	    maxband = ahivr (Memr[band], nwave)

	    if (maxband != 0.0) {
		do iwave = 1, nwave {
		    frac = Memr[band+iwave-1] / maxband
		    Memr[band+iwave-1] = minval * (1.0 - frac) + maxval * frac
		}
	    }

	    band = band + nwave
	}

end

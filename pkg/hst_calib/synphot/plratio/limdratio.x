include "dratio.h"

#* HISTORY *
#* B.Simon	12-Jul-94	original

# LIMDRATIO -- Compute plot limits for spectrum from flux limits

procedure limdratio (dratio, inlimit, outlimit)

pointer	dratio		# u: spectrum descriptor
real	inlimit[4]	# i: input plot limits
real	outlimit[4]	# u: calculated plot limits
#--
int	nwave, iflux
pointer	wave, flux

begin
	if (dratio == NULL)
	    return

	# Compute limits to spectrums

	wave = RAT_WAVE(dratio)
	nwave = RAT_NWAVE(dratio)
	flux = RAT_FLUX(dratio)

	do iflux = 1, RAT_NFLUX(dratio) {
	    call setlimit (nwave, Memr[wave], Memr[flux], inlimit, outlimit)
	    flux = flux + nwave
	}

end

include "dband.h"
include "../limit.h"

#* HISTORY *
#* B.Simon	08-Jun-94	original

# LIMDBAND -- Compute plot limits for bandpasses

procedure limdband (dband, inlimit, outlimit)

pointer	dband		# u: spectrum descriptor
real	inlimit[4]	# i: input plot limits
real	outlimit[4]	# u: calculated plot limits
#--
int	nwave, iband
pointer	wave, band

begin
	if (dband == NULL)
	    return

	# Compute limits to bandpasses

	wave = BND_WAVE(dband)
	nwave = BND_NWAVE(dband)
	band = BND_BAND(dband)

	do iband = 1, BND_NBAND(dband) {
	    call setlimit (nwave, Memr[wave], Memr[band], inlimit, outlimit)
	    band = band + nwave
	}

end

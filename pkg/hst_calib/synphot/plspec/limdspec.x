include "dspec.h"

#* HISTORY *
#* B.Simon	08-Jun-94	original

# LIMDSPEC -- Compute plot limits for spectrum from flux limits

procedure limdspec (dspec, form, inlimit, outlimit)

pointer	dspec		# u: spectrum descriptor
char	form[ARB]	# i: Output flux units
real	inlimit[4]	# i: input plot limits
real	outlimit[4]	# u: calculated plot limits
#--
int	nwave, iflux
pointer	wave, flux

begin
	if (dspec == NULL)
	    return

	# Compute limits to spectrums

	wave = SPC_WAVE(dspec)
	nwave = SPC_NWAVE(dspec)
	flux = SPC_FLUX(dspec)

	do iflux = 1, SPC_NFLUX(dspec) {
	    call fluxlimit (nwave, Memr[wave], Memr[flux], form, 
			    inlimit, outlimit)

	    flux = flux + nwave
	}

end

include "dsphot.h"

#* HISTORY *
#* B.Simon	13-Jun-94	original

# LIMDSPHOT -- Compute plot limits for spectrum from flux limits

procedure limdsphot (dsphot, inlimit, outlimit)

pointer	dsphot		# u: spectrum descriptor
real	inlimit[4]	# i: input plot limits
real	outlimit[4]	# u: calculated plot limits
#--
int	nwave, iflux
pointer	sp, wave, flux, err, fluxp, fluxm

begin
	if (dsphot == NULL)
	    return

	# Set pointers to arrays within structure

	wave = SPT_WAVE(dsphot)
	nwave = SPT_NWAVE(dsphot)
	flux = SPT_FLUX(dsphot)
	err = SPT_ERROR(dsphot)

	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (fluxp, nwave, TY_REAL)
	call salloc (fluxm, nwave, TY_REAL)

	# Compute limits to spectrums

	do iflux = 1, SPT_NFLUX(dsphot) {
	    if (IS_INDEFR (Memr[err])) {
		call setlimit (nwave, Memr[wave], Memr[flux], 
			      inlimit, outlimit)

	    } else {
		call aaddr (Memr[flux], Memr[err], Memr[fluxp], nwave)
		call asubr (Memr[flux], Memr[err], Memr[fluxm], nwave)

		call setlimit (nwave, Memr[wave], Memr[fluxp], 
			       inlimit, outlimit)
		call setlimit (nwave, Memr[wave], Memr[fluxm], 
			       inlimit, outlimit)
	    }

	    flux = flux + nwave
	    err = err + nwave
	}

	call sfree (sp)
end

include	"simtwo.h"

#* HISTORY *
#* B.Simon	19-Jun-95	original

# EARTHSHINE -- Compute background due to earth light

procedure earthshine (earthtab, eshine, exptime, apscale, nwave, wave, 
		      thruput, espec, ecount)

char	earthtab[ARB]	# i: earthlight spectrum
real	eshine		# i: fraction of maximum earthshine (0. - 1.)
real	exptime		# i: exposure time
double	apscale		# i: aperture scale
int	nwave		# i: length of wavelength and thruput tables
real	wave[ARB]	# i: wavelength array
real	thruput[ARB]	# i: observation mode thruput
real	espec[ARB]	# o: earthlight spectrum
real	ecount		# o: counts / pix from earth light
#--
pointer	sp, flux
real	pixlen, area, factor

real	effstim()

begin
	# Check for null file

	if (earthtab[1] == EOS) {
	    ecount = 0.0
	    return
	}

	call smark (sp)
	call salloc (flux, nwave, TY_REAL)

	# Compute pixel area in square arcseconds

	pixlen = 1.0
	call angtodegr (BACK_UNITS, pixlen)
	pixlen = apscale / pixlen
	area = pixlen * pixlen

	# Read in earthlight spectrum. Flux is assumed to be per sq sec.

	call rdspec (earthtab, nwave, wave, Memr[flux])

	# Multiply by factor

	factor = area * exptime * eshine
	call amulkr (Memr[flux], factor, Memr[flux], nwave)

	# Compute counts from observation mode throughput and flux

	ecount = effstim (nwave, wave, thruput, Memr[flux], OBJ_UNITS)

	# Compute earthlight spectrum

	call amulr (thruput, Memr[flux], espec, nwave)

	call sfree (sp)
end

include	"simtwo.h"

#* HISTORY *
#* B.Simon	19-Jun-95	original

# THERMAL -- Compute background due to earth light

procedure thermal (thermtab, exptime, apscale, nwave, wave, 
		      thruput, tcount)

char	thermtab[ARB]	# i: thermal spectrum
real	exptime		# i: exposure time
double	apscale		# i: aperture scale
int	nwave		# i: length of wavelength and thruput tables
real	wave[ARB]	# i: wavelength array
real	thruput[ARB]	# i: observation mode thruput
real	tcount		# o: counts / pix from earth light
#--
pointer	sp, flux
real	pixlen, area, factor

real	effstim()

begin
	# Check for null file

	if (thermtab[1] == EOS) {
	    tcount = 0.0
	    return
	}

	call smark (sp)
	call salloc (flux, nwave, TY_REAL)

	# Compute pixel area in square arcseconds

	pixlen = 1.0
	call angtodegr (BACK_UNITS, pixlen)
	pixlen = apscale / pixlen
	area = pixlen * pixlen

	# Read in thermal spectrum. Flux is assumed to be per sq sec.

	call rdspec (thermtab, nwave, wave, Memr[flux])

	# Multiply by factor

	factor = area * exptime
	call amulkr (Memr[flux], factor, Memr[flux], nwave)

	# Compute counts from observation mode throughput and flux

	tcount = effstim (nwave, wave, thruput, Memr[flux], OBJ_UNITS)

	call sfree (sp)
end

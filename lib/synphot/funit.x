include	"libsynphot.h"

# FUNIT -- Compute the stimulus needed to produce a unit response from the HST

#  This procedure computes the stimulus that is needed to produce a unit 
#  response of one count per second. This is a measure of the inverse 
#  sensitivity of the passband.
#
#  The unit stimulus is given in units of FLAM (erg s^-1 cm^-2 A^-1) by
#  the formula:
#
#   U_lam = hc / [ A * INT(thru * lam * dlam) ]
 
real procedure funit (area, nwave, wave, thruput)

real	area		# i: telescope area
int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavelength array
real	thruput[ARB]	# i: throughput array
#--
pointer	sp, fun
real	unit

real	syntegral()

begin
	call smark (sp)
	call salloc (fun, nwave, TY_REAL)

	call aabsr (thruput, Memr[fun], nwave)
	call amulr (wave, Memr[fun], Memr[fun], nwave)

	unit = syntegral (nwave, wave, Memr[fun])

	if (unit <= 0.0) {
	    unit = 0.0
	} else {
	    unit = (H * C) / (area * unit)
	}

	call sfree (sp)
	return (unit)
end

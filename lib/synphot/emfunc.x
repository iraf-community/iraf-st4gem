include "libsynphot.h"

# EMFUNC -- Compute a gaussian emission line spectrum 

procedure emfunc (mean, fwhm, flux, units, nwave, wave, spec)

real	mean		# i: mean wavelength of emission line
real	fwhm		# i: full width at half maximum
real	flux		# i: total flux in emission line
char	units[ARB]	# i: flux units
int	nwave		# i: length of wavelength and passband arrays
real	wave[ARB]	# i: wavelengths at which spectrum is computed
real	spec[ARB]	# i: emission line spectrum
#--
int	utype, iwave, done
pointer	sp, form
real	sigma, factor, area

string	formlist FORMSTR

int	word_match(), anytophot()

begin
	# Allocate memory for temporary string

	call smark (sp)
	call salloc (form, SZ_FNAME, TY_CHAR)

	# Compute standard deviation and normalization factor

	sigma = fwhm / sqrt (8.0 * log (2.0))
	factor = flux / (sqrt (2.0 * PI) * sigma)

	# Check for count units and correct for hst area

	call strcpy (units, Memc[form], SZ_FNAME)
	call strfix (Memc[form])

	utype = word_match (Memc[form], formlist)

	switch (utype) {
	case 2:		# counts
	    call strcpy ("photlam", Memc[form], SZ_FNAME)
	    call get_hstarea (area)
	    factor = factor / area

	case 11:	# obmag
	    call strcpy ("photlam", Memc[form], SZ_FNAME)
	    call get_hstarea (area)
	    factor = 10 ** (-0.4 * factor) / area
	}

	# Compute emission line with gaussian profile

	do iwave = 1, nwave
	    spec[iwave] = factor * exp (-0.5 *((wave[iwave] - mean)/ sigma)**2)

	# Convert to photlam units

	done = anytophot (Memc[form], nwave, wave, spec)
	call sfree (sp)
end

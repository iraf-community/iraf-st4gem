include "libsynphot.h"

# RATE -- Compute the count rate produced by a spectrum in a passband

real procedure rate (nwave, wave, band, spec, units)

int	nwave		# i: length of passband and spectral arrays
real	wave[ARB]	# i: wavelength at which spectrum is computed
real	band[ARB]	# i: instrument bandpass
real	spec[ARB]	# i: spectrum
char	units[ARB]	# i: units of spectrum
#--
int	utype, done
pointer	sp, prod, form
real	result

string	badunits "Unknown flux units"
string	formlist FORMSTR

int	is_magunit(), word_match(), anytophot(), phottoany()
real	sumfilt(), asumr()

errchk	anytophot, phottoany, synphoterr

begin
	# Allocate memory for temporary arrays

        call smark (sp)
	call salloc (prod, nwave, TY_REAL)
        call salloc (form, SZ_FNAME, TY_CHAR)

	# Compute product of bandpass and spectrum

	if (is_magunit (units) == NO) {
	    call amulr (band, spec, Memr[prod], nwave)

	} else {
	    call amovr (spec, Memr[prod], nwave)
	    done = anytophot (units, nwave, wave, Memr[prod])

	    call amulr (band, Memr[prod], Memr[prod], nwave)
	}

	# Compute count rate by integrating over product

        call strcpy (units, Memc[form], SZ_FNAME)
        call strfix (Memc[form])

	utype = word_match (Memc[form], formlist)

	switch (utype) {
	case 1: # photlam
	    result = sumfilt (nwave, wave, 0, Memr[prod])

	case 2: # counts
	    result = asumr (Memr[prod], nwave)

	case 3: # flam
	    result = sumfilt (nwave, wave, 1, Memr[prod])
	    result = result / (H * C)

	case 4: # fnu
	    result = sumfilt (nwave, wave, -1, Memr[prod])
	    result = result / H

	case 5: # photnu
	    result = sumfilt (nwave, wave, -2, Memr[prod])
	    result = result * C

	case 6: # jy
	    result = sumfilt (nwave, wave, -1, Memr[prod])
	    result = result * 1.0e-23 / H

	case 7: # mjy
	    result = sumfilt (nwave, wave, -1, Memr[prod])
	    result = result * 1.0e-26 / H

	case 8: # abmag
	    done = phottoany ("fnu", nwave, wave, Memr[prod])
	    result = sumfilt (nwave, wave, -1, Memr[prod])
	    result = result / H

	case 9: # stmag
	    done = phottoany ("flam", nwave, wave, Memr[prod])
	    result = sumfilt (nwave, wave, 1, Memr[prod])
	    result = result / (H * C)

	case 10: # vegamag
	    done = phottoany ("flam", nwave, wave, Memr[prod])
	    result = sumfilt (nwave, wave, 1, Memr[prod])
	    result = result / (H * C)

	case 11: # obmag
	    done = phottoany ("counts", nwave, wave, Memr[prod])
	    result = asumr (Memr[prod], nwave)

	default: # unknown units
	    call synphoterr (badunits, units)
	}

	call sfree (sp)
	return (result)
end

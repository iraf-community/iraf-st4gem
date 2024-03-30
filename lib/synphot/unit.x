include "libsynphot.h"

# UNIT -- Compute the stimulus needed for unit response per cm^2 sec

real procedure unit (nwave, wave, band, units)

int	nwave		# i: length of passband and spectral arrays
real	wave[ARB]	# i: wavelength at which spectrum is computed
real	band[ARB]	# i: instrument bandpass
char	units[ARB]	# i: units of spectrum
#--
int	utype
pointer	sp, form,vflux
real	sum, resp

string	badunits "Unknown flux units"
string	formlist FORMSTR

int	word_match()
real	sumfilt(), asumr(), rate()

errchk	rdstospec, rate, synphoterr

begin
	# Allocate memory for temporary arrays

        call smark (sp)
        call salloc (form, SZ_FNAME, TY_CHAR)

	# Initialize output

	resp = 0.0

	# Compute unit response by integrating over bandpass

        call strcpy (units, Memc[form], SZ_FNAME)
        call strfix (Memc[form])

	utype = word_match (Memc[form], formlist)

	switch (utype) {
	case 1: # photlam
	    sum = sumfilt (nwave, wave, 0, band)
	    if (sum > 0.0)
		resp = 1.0 / sum

	case 2: # counts
	    sum = asumr (band, nwave)
	    if (sum > 0.0)
		resp = 1.0 / sum

	case 3: # flam
	    sum = sumfilt (nwave, wave, 1, band)
	    sum = sum / (H * C)
	    if (sum > 0.0)
		resp = 1.0 / sum

	case 4: # fnu
	    sum = sumfilt (nwave, wave, -1, band)
	    sum = sum / H
	    if (sum > 0.0)
		resp = 1.0 / sum

	case 5: # photnu
	    sum = sumfilt (nwave, wave, -2, band)
	    sum = sum * C
	    if (sum > 0.0)
		resp = 1.0 / sum

	case 6: # jy
	    sum = sumfilt (nwave, wave, -1, band)
	    sum = sum * 1.0e-23 / H
	    if (sum > 0.0)
		resp = 1.0 / sum

	case 7: # mjy
	    sum = sumfilt (nwave, wave, -1, band)
	    sum = sum * 1.0e-26 / H
	    if (sum > 0.0)
		resp = 1.0 / sum

	case 8: # abmag
	    sum = sumfilt (nwave, wave, -1, band)
	    sum = sum / H
	    if (sum > 0.0)
		resp = 2.5 * alog10 (sum) + ABZERO

	case 9: # stmag
	    sum = sumfilt (nwave, wave, 1, band)
	    sum = sum / (H * C)
	    if (sum > 0.0)
		resp = 2.5 * alog10 (sum) + STZERO

	case 10: # vegamag
	    call salloc (vflux, nwave, TY_REAL)
	    call rdstospec (VEGA, nwave, wave, Memr[vflux])
	    sum = rate (nwave, wave, band, Memr[vflux], "photlam")
	    if (sum > 0.0)
		resp = 2.5 * alog10 (sum)

	case 11: # obmag
	    sum = asumr (band, nwave)
	    if (sum > 0.0)
		resp = 2.5 * alog10 (sum)

	default: # unknown units
	    call synphoterr (badunits, units)
	}

	call sfree (sp)
	return (resp)
end

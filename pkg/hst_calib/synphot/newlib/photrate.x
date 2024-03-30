#* HISTORY *
#* B.Simon	15-Jun-94	original

# PHOTRATE -- Compute the photometric rate for a unit flux

real procedure photrate (form, wave, nwave, band1, band2, nband)

char	form[ARB]	# i: spectral form
real	wave[ARB]	# i: wavelength array
int	nwave		# i: number of wavelengths
real	band1[ARB]	# i: first bandpass
real	band2[ARB]	# i: second bandpass (or zero)
int	nband		# i: number of bandpasses
#--
pointer	sp, flux, prod
real	factor, area, base, factor2

string	badform    "Unrecognized form"
string	badfactor  "Photometric rate is zero"

int	is_count(), is_magunit(), anytophot()
real	syntegral()

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (flux, nwave, TY_REAL)
	call salloc (prod, nwave, TY_REAL)

	if (is_count (form) == YES) {
	    # The photometric rate for counts is the inverse telescope area

	    call get_hstarea (area)
	    factor = 1.0 / area

	} else {
	    # Convert unix flux in photlam

	    if (is_magunit (form) == NO) {
		base = 1.0
	    } else {
		base = 0.0
	    }

	    call amovkr (base, Memr[flux], nwave)
	    if (anytophot (form, nwave, wave, Memr[flux]) == NO)
		call printerr_str (badform, form)

	    # Integrate over bandpass to get photometric rate

	    call amulr (band1, Memr[flux], Memr[prod], nwave)
	    factor = syntegral (nwave, wave, Memr[prod])

	    if (nband > 1) {
		call amulr (band2, Memr[flux], Memr[prod], nwave)
		factor2 = syntegral (nwave, wave, Memr[prod])
		if (factor2 > 0) {
		    factor = factor / factor2
		} else {
		    call printerr_real (badfactor, factor2)
		}
	    }
	}

	# Return photometric rate for a unit stimulus

	call sfree (sp)
	return (factor)
end

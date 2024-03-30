#* HISTORY *
#* B.Simon	17-May-94	Adapted from renormfunc
#* B.Simon	27-Apr-95	Changed test for zero stimulus

# EFFSTIM -- Compute the effective stimulus of passband and spectrum

real procedure effstim (nwave, wave, filt, spec, form)

int	nwave		# i: length of wavelength and spectral arrays
real	wave[ARB]	# i: wavelengths at which spectrum is computed
real	filt[ARB]	# i: bandpass to normalize over
real	spec[ARB]	# i: spectrum to be normalized
char	form[ARB]	# i: units of constant value
#--
int	done
pointer	sp, spec2
real	cntrate, resp, result

string	blankform  "Form is blank"

int	phottoany(), is_magunit(), is_count()
real	rate(), unit()

begin
	call smark (sp)
	call salloc (spec2, nwave, TY_REAL)

	# Convert spectrum to renormalization form

	call amovr (spec, Memr[spec2], nwave)
	done = phottoany (form, nwave, wave, Memr[spec2])

	if (done == NO)
	    call synphoterr (blankform, form)

	# Compute count rate and unit response for this passband and spectrum

	cntrate = rate (nwave, wave, filt, Memr[spec2], form)
	resp = unit (nwave, wave, filt, form)

	# If the units are in counts, the effective stimulus is the count rate
	# Otherwise, the effective stimulus is the product of the count rate
	# and the unit response. 

	if (is_magunit (form) == NO) {
	    if (is_count (form) == NO) {
		if (cntrate <= 0.0 || resp <= 0.0) {
		    result = 0.0
		} else {
		    result = cntrate * resp
		}

	    } else {
		if (cntrate <= 0.0) {
		    result = 0.0
		} else {
		    result = cntrate
		}
	    }

	} else {
	    if (is_count (form) == NO) {
		if (cntrate <= 0.0) {
		    result = 100.0
		} else {
		    result = resp - 2.5 * alog10 (cntrate)
		}

	    } else {
		if (cntrate <= 0.0) {
		    result = 100.0
		} else {
		    result = -2.5 * alog10 (cntrate)
		}
	    }
	}

	call sfree (sp)
	return (result)
end

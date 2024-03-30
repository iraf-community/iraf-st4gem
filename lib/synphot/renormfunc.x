#* HISTORY *
#* B.Simon	27-Apr-95	Changed test for zero stimulus
#* B.Simon	26-Jan-96	Removed iferr around renormfunc()

# RENORMFUNC -- Renormalize the flux of a spectrum

procedure renormfunc (inspec, band, const, units, nwave, wave, spec)

real	inspec[ARB]	# i: spectrum to be normalized
real	band[ARB]	# i: bandpass to normalize over
real	const		# i: value to normalize to
char	units[ARB]	# i: units of constant value
int	nwave		# i: length of wavelength and spectral arrays
real	wave[ARB]	# i: wavelengths at which spectrum is computed
real	spec[ARB]	# o: normalized spectrum, in photlam
#--
int	done
pointer	sp, outspec
real	cntrate, resp, effstim, factor

string	badunits  "Unknown flux units for function renorm"
string	badrenorm "Cannot renormalize spectrum, flux is zero in passband"

int	phottoany(), anytophot(), is_magunit(), is_count()
real	rate(), unit()

errchk	rate, unit, anytophot, synphoterr

begin
	call smark (sp)
	call salloc (outspec, nwave, TY_REAL)

	# Convert spectrum to renormalization units

	call amovr (inspec, Memr[outspec], nwave)

	done = phottoany (units, nwave, wave, Memr[outspec])
	if (done == NO)
	    call synphoterr (badunits, units)

	# Compute count rate and unit response for this passband and spectrum

	cntrate = rate (nwave, wave, band, Memr[outspec], units)
	resp = unit (nwave, wave, band, units)

	# Compute renormalization factor

	if (cntrate <= 0.0)
	    call synphoterr (badrenorm, "rn")

	if (is_magunit (units) == NO) {
	    if (is_count(units) == NO) {
		if (resp <= 0.0)
		    call synphoterr (badrenorm, "rn")

		effstim = cntrate * resp
	    } else {
		effstim = cntrate
	    }

	    factor = const / effstim
	    call amulkr (Memr[outspec], factor, spec, nwave)

	} else {
	    if (is_count(units) == NO) {
		effstim = resp - 2.5 * alog10 (cntrate)
	    } else {
		effstim = -2.5 * alog10 (cntrate)
	    }

	    factor = const - effstim
	    call aaddkr (Memr[outspec], factor, spec, nwave)
	}

	# Convert spectrum back to photlam units

	done = anytophot (units, nwave, wave, spec)
	call sfree (sp)
end

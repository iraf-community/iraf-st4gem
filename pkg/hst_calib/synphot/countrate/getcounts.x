# GETCOUNTS -- Compute counts at reference wavelength

#* HISTORY *
#* B.Simon	17-Mar-93	original
#* B.Simon	15-Apr-94	conversion of spectrum to counts added
#* B.Simon	27-May-94	removed conversion of units to counts

real procedure getcounts (form, refwave, nwave, wave, flux)

char	form[ARB]	# i: flux units
real	refwave		# i: reference wavelength
int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavelength array
real	flux[ARB]	# i: spectrum as a function of wavelength
#--
int	iw, jw
real	p, count_ref

begin
	# Interpolate for counts at reference wavelength
	# Divide by wavelength range to get counts per angstrom

	if (IS_INDEFR (refwave)) {
	    count_ref = INDEFR

	} else {

	    # Find wavelengths that bracket reference wavelength

	    jw = 0
	    do iw = 1, nwave {
		if (wave[iw] > refwave) {
		    jw = iw
		    break
		}
	    }

	    # Interpolate for value at wavelength

	    if (jw < 2) {
		count_ref = INDEFR

	    } else {
		p = (refwave - wave[jw-1]) / (wave[jw] - wave[jw-1])
		count_ref = flux[jw-1] * (1.0 - p) + flux[jw] * p
	    }
	}

	return (count_ref)
end

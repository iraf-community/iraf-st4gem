#* HISTORY *
#* B.Simon	08-Jun-94	adopted from photstim

# EFFSTIM2 -- Compute the effective stimulus for calcphot

real procedure effstim2 (form, nwave, wave, nband, band1, band2, flux)

char	form[ARB]	# i: Form or type of result
int	nwave		# i: length of wavelength set
real	wave[ARB]	# i: wavelength set
int	nband		# i: Number of bandpasses
real	band1[ARB]	# i: First bandpass 
real	band2[ARB]	# i: Second bandpass (or zero)
real	flux[ARB]	# i: Spectral flux
#--
real	result, result1, result2

string	nostim   "Effective stimulus is zero or negative"

int	is_magunit()
real	effstim()

begin
	result1 = effstim (nwave, wave, band1, flux, form)

	# If nband > 1, compute effective stimulus 
	# for each passband and combine the results

	if (nband == 1) {
	    result = result1

	} else {
	    result2 = effstim (nwave, wave, band2, flux, form)

	    if (is_magunit (form) == YES) {
		result = result1 - result2
	    } else {
		if (result2 > 0.0) {
		    result = result1 / result2
		} else {
		    call printerr_real (nostim, result2)
		}
	    }
	}

	return (result)
end


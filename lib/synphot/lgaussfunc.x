# LGAUSSFUNC -- Compute a gaussian passband

procedure lgaussfunc (mean, fwhm, nwave, wave, band)

real	mean		# i: mean of gaussian function
real	fwhm		# i: full width at half maximum
int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelengths at which bandpass is computed
real	band[ARB]	# o: instrument bandpass
#--
int	iwave
real	sigma, arg

string	badwidth  "Width must be a positive number"

begin
	if (fwhm <= 0.0)
	    call synphoterr (badwidth, "lgauss")

	arg = (fwhm * fwhm) / (4.0 * mean * mean)
	if (arg <= 0.1) {
	    sigma = arg /  sqrt (8.0 * log(2.0))
	} else {
	    arg = 1.0 - arg
	    sigma = - log (arg) / sqrt (8.0 * log(2.0))
	}

	do iwave = 1, nwave
	    band[iwave] = exp (-0.5 *(log(wave[iwave] / mean)/ sigma)**2)

end

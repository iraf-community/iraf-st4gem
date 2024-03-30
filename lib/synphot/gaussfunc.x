# GAUSSFUNC -- Compute a gaussian passband

procedure gaussfunc (mean, fwhm, nwave, wave, band)

real	mean		# i: mean of gaussian function
real	fwhm		# i: full width at half maximum
int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelengths at which bandpass is computed
real	band[ARB]	# o: instrument bandpass
#--
int	iwave
real	sigma

string	badwidth  "Width must be a positive number"

begin
	if (fwhm <= 0.0)
	    call synphoterr (badwidth, "gauss")

	sigma = fwhm / sqrt (8.0 * log(2.0))

	do iwave = 1, nwave
	    band[iwave] = exp (-0.5 *((wave[iwave] - mean)/ sigma)**2)

end

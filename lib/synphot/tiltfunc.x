# TILTFUNC -- Multiply a passband by a Legendre polynomial

procedure tiltfunc (pass, npoly, poly, nwave, wave, band)

real	pass[ARB]	# i: input passband
int	npoly		# i: number of polynomial coefficients
real	poly[ARB]	# i: polynomial coefficients
int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelengths at which bandpass is computed
real	band[ARB]	# o: instrument bandpass
#--
real	mean, fwhm

real	avglam(), fwhmlam()
errchk	polyfunc

begin
	# Compute mean wavelength and width from input wavelength

	mean = avglam (nwave, wave, pass)
	fwhm = fwhmlam (nwave, wave, pass)

	# Calculate legendre polynomial and multiply it 
	# by input passband to generate output passband

	call polyfunc (mean, fwhm, npoly, poly, nwave, wave, band)
	call amulr (pass, band, band, nwave)

end

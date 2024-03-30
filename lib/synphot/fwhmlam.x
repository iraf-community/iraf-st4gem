# FWHMLAM -- Compute the full width half maximum wavelength

real procedure fwhmlam (nwave, wave, thruput)

int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavelength array
real	thruput[ARB]	# i: throughput array
#--
real	fwhm
real	rmslam()

begin
	fwhm = sqrt (8.0 * log (2.0)) * rmslam (nwave, wave, thruput)
	return (fwhm)
end

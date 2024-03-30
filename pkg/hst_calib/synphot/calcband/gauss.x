define	FACTOR	sqrt( 8.*log(2.) )

# GAUSS -- Create a gaussian passband normalized to maximum value of 1

procedure gauss(script, iw, nwave, wave, band)

char	script[ARB]	# i: Command script
int	iw		# io: Position in script
int	nwave		# i: Number of wavelengths in wave
real	wave[ARB]	# i: Wavelength array
real	band[ARB]	# o: Passband

int	nchar, ic
real	mu, sigma, fwhm
real	ctor()

begin

	#Get average

	nchar = ctor( script, iw, mu)

	#Get fwhm:  full width at half maximum

	nchar = ctor( script, iw, fwhm)
	sigma = fwhm/FACTOR

	# Calculate gaussian

	do ic = 1, nwave
	   band[ic] = exp( - 0.5 * ( (wave[ic] - mu)/sigma )**2 )

end

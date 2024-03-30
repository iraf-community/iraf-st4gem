define	FACTOR	sqrt( 8.*log(2.) )

# LGAUSS -- Create a gaussian passband in log wavelength, normalized
#           to maximum value of 1

procedure lgauss(script, iw, nwave, wave, band)

char	script[ARB]	# i: Command script
int	iw		# io: Position in script
int	nwave		# i: Number of wavelengths in wave
real	wave[ARB]	# i: Wavelength array
real	band[ARB]	# o: Passband

int	nchar, ic
real	mu, sigma, fwhm, arg
int	ctor()

begin

	#Get average

	nchar = ctor( script, iw, mu)

	#Get fwhm:  full width at half maximum

	nchar = ctor( script, iw, fwhm)

	arg = fwhm*fwhm/(4.*mu*mu)
	if (arg <= 0.1)
	   sigma = arg/FACTOR

	else {
	   arg = 1. - arg
	   sigma = - log(arg)/FACTOR
	}

	# Calculate gaussian

	do ic = 1, nwave
	   band[ic] = exp( - 0.5 * ( log(wave[ic]/mu)/sigma )**2 )

end

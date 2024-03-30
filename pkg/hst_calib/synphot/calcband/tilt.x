define	FACTOR		sqrt(8.*log(2.))
define	MAXORDER	10
define	MAXARGVAL	69.

# TILT -- Routine to multiply a passband by a polynomial
#         in u where u = (x - x0)/sigma.  x0 is average wavelength of 
#         passband sigma is half width at half maximum.

procedure tilt( script, nwave, wave, mu, fwhm, ip, band )

char	script[ARB]	# i: Command script
int	nwave		# i: Number of wavelength in wave
real	wave[ARB]	# i: Wavelength array
real	mu		# i: Offset for polynomial
real	fwhm		# i: Scaling for polynomial
int	ip		# io: Position in script
real	band[ARB]	# io: Passband array

int	ic, iw, ncoeff
real	val, poly, u, coeff[MAXORDER], sigma, funct
int	ctor()
real	legendre()

begin

	# Get polynomial coeffecients from script.  Stop at non-number
	# Maximum order of 10

	for (ic = 1; ctor(script, ip, val) > 0 && ic <= MAXORDER+1; ic = ic+1)
	   coeff[ic] = val

	ncoeff = ic - 1
	sigma = fwhm/FACTOR

	do iw = 1,nwave {

	   poly = 0.
	   u = (wave[iw] - mu)/sigma

	   # Calculate the value of the polynomial at the current wavelength
	   do ic = 1, ncoeff
	      poly = poly + coeff[ic] * legendre( u, ic )

	   if ( poly >= 0 )
	      funct = 1. +  poly
	   else
	      funct = exp( max( poly, -MAXARGVAL) )

	   if ( !IS_INDEFR (band[iw]) )
	      band[iw] = band[iw] * funct
	   else
	      band[iw] = INDEFR

	}

end

# EBMVLFUNC -- Evaluation of Howarth's LMC extinction curve
#
# Values of LMC extinction are from Howarth's 1983 paper MNRAS, 203, 301.
#
# Extinction:
# x = 1/lambda(microns) = 10000/lambda(angstroms)
# X(x) = A(lambda)/E(B-V)
# R = A(V)/E(B-V) = 3.1
#
# Infrared
#     x <= 1.83
#     X(x) = [(1.86 - 0.48*x)*x - 0.1]*x
#
# Visible (Optical)
#     1.83 <= x <= 2.75
#     X(x) = R + 2.04*(x - 1.83) + 0.094*(x - 1.83)**2
#
# Ultra Violet
#     2.75 <= x <= 9.0
#     X(x) = R - 0.236 + 0.462*x + 0.105*x*x + 0.454/[(x - 4.557)**2 + 0.293]

procedure ebmvlfunc (extval, nwave, wave, band)

real	extval		# i: extinction value
int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelength set output is produced on
real	band[ARB]	# o: output bandpass
#--
int	iwave
real	x, extl

begin
	
	do iwave = 1, nwave {
	    # Convert wavelength in angstroms to 1/microns
	    x = 10000.0 / wave[iwave]

	    # Infrared - extend optical results linearly to 0 at 1/lam = 0
	    if ( x <= 1.83)
		extl = ((1.86 - 0.48 * x) * x - 0.1) * x

	    else if ( x <= 2.75 ) { 
		extl = 3.1 + 2.04 * (x - 1.83) + 0.094 * 
		       (x - 1.83) * (x - 1.83)
		
	    # Continue out to lambda = 912 A
	    } else { 
		x = min (x, 10.96)
		extl = 3.1 - 0.236 + 0.462 * x + 0.105 * x * x +
		       0.454 / ((x - 4.557) * (x - 4.557) + 0.293)
	    }

	    extl = extl * extval
	    band[iwave] = 10.0 ** (-0.4 * extl)
	}
end

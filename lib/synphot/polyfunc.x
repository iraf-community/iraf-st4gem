define	MAXEXPON	69.0

# POLYFUNC -- Compute a passband from a Legendre polynomial

procedure polyfunc (mean, fwhm, npoly, poly, nwave, wave, band)

real	mean		# i: mean wavelength
real	fwhm		# i: full width at half maximum
int	npoly		# i: number of polynomial coefficients
real	poly[ARB]	# i: polynomial coefficients
int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelengths at which bandpass is computed
real	band[ARB]	# o: instrument bandpass
#--
int	iwave, ipoly
real	sigma, x, u

real	legfunc()

begin
	# Convert full width, half maximum to standard deviation

	sigma = fwhm / sqrt (8.0 * log(2.0))

	do iwave = 1, nwave {
	    # Compute legendre polynomials

	    u = (wave[iwave] - mean) / sigma

	    # Multiply by input coefficients to get value

	    x = 0.0
	    do ipoly = 1, npoly 
		x = x + poly[ipoly] * legfunc (u, ipoly)

	    # Convert function to final form

	    if (x >= 0.0) {
		band[iwave] = x + 1.0
	    } else {
		x = max (x, - MAXEXPON)
		band[iwave]= exp (x)
	    }
	}

end

# LEGFUNC -- Compute Legendre polynomial

real procedure legfunc (u, npoly) 

real	u		# i: independent variable
int	npoly		# i: degree of polynomial
#--
int	ipoly
pointer	sp, leg

begin
	call smark (sp)
	call salloc (leg, npoly+2, TY_REAL)

	Memr[leg] = 1.0
	Memr[leg+1] = u

	do ipoly = 2, npoly {
	    Memr[leg+ipoly] = (2 * ipoly - 1) * u * Memr[leg+ipoly-1] -
			      (npoly - 1) * Memr[leg+ipoly-2]

	    Memr[leg+ipoly] = Memr[leg+ipoly] / ipoly
	}

	call sfree (sp)
	return (Memr[leg+npoly])
end

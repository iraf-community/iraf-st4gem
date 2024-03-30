#* HISTORY *
#* B.Simon	08-Mar-93	original

# SPECINTRP -- Interpolate a spectrum onto a new grid

procedure specintrp (iwave, ispec, ilen, owave, ospec, olen)

real	iwave[ARB]	# i: input wavelength grid
real	ispec[ARB]	# i: input flux
int	ilen		# i: length of input spectrum
real	owave[ARB]	# i: output wavelength grid
real	ospec[ARB]	# o: output flux
int	olen		# i: length of output spectrum
#--
int	ipix, opix
real	ilo, ihi, olo, ohi, dwave, frac

begin

	# Compute the value of the output spectrum at each point 
	# by integrating over the flux of the input spectrum

	ipix = 1
	opix = 1

	while (opix <= olen) {

	    # Set the limits of integration

	    if (opix == 1) {
		olo = owave[1] - 0.5 * (owave[2] - owave[1])
		ohi = owave[1] + 0.5 * (owave[2] - owave[1])
	    } else {
		olo = ohi
		ohi = owave[opix] + 0.5 * (owave[opix] - owave[opix-1])
	    }


	    # Sum the input counts over the limits of integration

	    ospec[opix] = 0.0

	    while (ipix <= ilen) {
		if (! IS_INDEFR(ispec[ipix])) {

		    # Compute the limits of the input pixel

		    if (ipix == 1) {
			dwave = iwave[2] - iwave[1]
		    } else {
			dwave = iwave[ipix] - iwave[ipix-1]
		    }

		    ilo = iwave[ipix] - 0.5 * dwave
		    ihi = iwave[ipix] + 0.5 * dwave

		    # Compute fractional part of input pixel within
		    # limits of integration

		    frac = 1.0
		    if (ilo < olo)
			frac = frac - (olo - ilo) / dwave

		    if (ihi > ohi)
			frac = frac - (ihi - ohi) / dwave

		    # Add pixel or fractional part to sum

		    if (frac > 0.0)
			ospec[opix] = ospec[opix] + frac * ispec[ipix]

		    # Sum is complete if next pixel is not inside
		    # the limits of integration

		    if (ihi >= ohi)
			break
		}

		ipix = ipix + 1
	    }

	    opix = opix + 1
	}

end

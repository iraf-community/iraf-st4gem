# Copyright restrictions apply - see st4gem$copyright.stsdas 
# 
#* HISTORY *
#* B.Simon	08-Mar-93	original
#* B.Simon	20-Aug-96	modified to use photlam units and coef array
#* B.Simon	17-Dec-96	modified to use trapezoid integration
#* B.Simon	22-Apr-97	compute offset for coefficient array

# REBINSPEC -- Rebin a spectrum onto a new wavelength grid

procedure rebinspec (nwave1, wave1, flux1, coef, nwave2, wave2, flux2)

int	nwave1		# i: length of input spectrum
real	wave1[ARB]	# i: input wavelength grid
real	flux1[ARB]	# i: input flux
char	coef[ARB]	# i: optional wavelength coefficients for output
int	nwave2		# i: length of output spectrum
real	wave2[ARB]	# i: output wavelength grid
real	flux2[ARB]	# o: output flux
#--
double	lo2, hi2, w[3]
int	iw1, iw2, jw, kw, done, mwave
real	hiflux, loflux, area, a, b, c

int	stepoffset(), setinform()
real	one_interp()

begin
	# Check for wavelength coeffcients
	if (coef[1] == '(') {
	    call calcstep (coef, a, b, c, mwave)

	    # Compute the offset from the start of the coefficients
	    # to the start of the wave2 array

	    kw = stepoffset (wave2[1], a, b, c)

	} else {
	    a = INDEFR
	    b = INDEFR
	    c = INDEFR
	}

	# Compute the value of the output spectrum at each point 
	# by integrating over the flux of the input spectrum

	iw1 = 1
	iw2 = 1

	while (iw2 <= nwave2) {

	    # Set the limits of integration

	    if (! IS_INDEFR (a)) {
		do jw = 1, 3
		    w[jw] = ((double(a) * (iw2 + jw + kw - 3)) + double(b)) * 
			    (iw2 + jw + kw - 3) + double(c)

		lo2 = 0.5 * (w[1] + w[2])
		hi2 = 0.5 * (w[2] + w[3])

	    } else {
		if (iw2 == 1) {
		    lo2 = wave2[1] - 0.5 * (wave2[2] - wave2[1])
		    hi2 = 0.5 * (wave2[2] + wave2[1])

		} else if (iw2 == nwave2) {
		    lo2 = hi2
		    hi2 = 0.5 * (wave2[nwave2] - wave2[nwave2-1]) +
			  wave2[nwave2]

		} else {
		    lo2 = hi2
		    hi2 = 0.5 * (wave2[iw2+1] + wave2[iw2])
		}

	    }

	    # Compute flux at endpoints

	    # This function updates the value of iw1 so that wave1[iw1] 
	    # is just greater than the wavelength passed as the second arg
	    # If the condition cannot be met, iw1 is set to nwave1+1

	    hiflux = one_interp (iw1, hi2, nwave1, wave1, flux1)
	    loflux = one_interp (iw1, lo2, nwave1, wave1, flux1)

	    # Sum the input counts over the limits of integration

	    if (iw1 > nwave1) {
		flux2[iw2] = 0.5 * (hi2 - lo2) * (hiflux + loflux)

	    } else if (hi2 < wave1[iw1]) {
		flux2[iw2] = 0.5 * (hi2 - lo2) * (hiflux + loflux)

	    } else {
		flux2[iw2] = 0.5 * (wave1[iw1] - lo2) * (flux1[iw1] + loflux)
		iw1 = iw1 + 1

		while (iw1 <= nwave1) {
		    if (wave1[iw1] >= hi2)
			break

		    flux2[iw2] = 0.5 * (wave1[iw1] - wave1[iw1-1]) *
				 (flux1[iw1] + flux1[iw1-1]) + flux2[iw2]
		    iw1 = iw1 + 1
		}

		flux2[iw2] = 0.5 * (hi2 - wave1[iw1-1]) * 
			     (hiflux + flux1[iw1-1]) + flux2[iw2]
	    }

	    iw2 = iw2 + 1
	}

	# Convert integrated flux to counts by multiplying by telescope area

	call gethstarea (area)
	call amulkr (flux2, area, flux2, nwave2)

	# Convert flux from counts to photlam units

	done = setinform ("counts", coef, nwave2, wave2, flux2)

end

#* HISTORY *
#* B.Simon	18-Jul-96	original
#* B.Simon	20-Aug-96	pass coef array instead of wavecat
#* B.Simon	22-Apr-97	compute offset for coefficient array

# A special routine for units conversion is called if the 
# units are counts and the wavelength set was specified by a set
# of coefficients. These routines allow more accurate computation
# of the width of each wavelength bin , which is very important
# for the highly sampled echelle modes.

# SETINFORM -- Convert the input flux into photlam units

procedure setinform (form, coef, nwave, wave, flux)

char	form[ARB]	# i: output form
char	coef[ARB]	# i: optional wavelength coefficients
int	nwave		# i: length of wavelength and flux arrays
real	wave[ARB]	# i: array of wavelengths
real	flux[ARB]	# u: array of fluxes
#--
int	done

int	is_count(), anytophot()

begin
	# Check for count units

	if (is_count (form) == NO) {
	    done = anytophot (form, nwave, wave, flux)

	} else {
	    # Check for wavelengths specified by coefficient set

	    if (coef[1] == '(') {
		call cnttophot (form, coef, nwave, wave, flux)
	    } else {
		done = anytophot (form, nwave, wave, flux)
	    }
	}

end

# SETOUTFORM -- Set the flux from photlam to the user specified output form

procedure setoutform (form, coef, nwave, wave, flux)

char	form[ARB]	# i: output form
char	coef[ARB]	# i: optional wavelength coefficients
int	nwave		# i: length of wavelength and flux arrays
real	wave[ARB]	# i: array of wavelengths
real	flux[ARB]	# u: array of fluxes
#--
int	done

int	is_count(), phottoany()

begin
	# Check for count units

	if (is_count (form) == NO) {
	    done = phottoany (form, nwave, wave, flux)

	} else {
	    # Check for wavelengths specified by coefficient set

	    if (coef[1] == '(') {
		call phottocnt (form, coef, nwave, wave, flux)
	    } else {
		done = phottoany (form, nwave, wave, flux)
	    }
	}

end

# CNTTOPHOT -- Special conversion routine for count to photlam units

procedure cnttophot (form, coef, nwave, wave, flux)

char	form[ARB]	# i: output form
char	coef[ARB]	# i: wavelength set coefficients
int	nwave		# i: length of wavelength and flux arrays
real	wave[ARB]	# i: array of wavelengths
real	flux[ARB]	# u: array of fluxes
#--
int	iwave, mwave, off
real	a, b, c, area, dwave

int	stepoffset(), is_magunit()

begin
	# The coefficients define a quadratic function used to
	# define the wavelength set. The wavelength spacing is
	# computed from the derivative of the quadratic

	call gethstarea (area)
	call calcstep (coef, a, b, c, mwave)
	off = stepoffset (wave[1], a, b, c)

	# If the form is not in magnitude units, it is COUNTS
	# If it is, it is OBMAG

	if (is_magunit (form) == NO) {
	    do iwave = 1, nwave {
		dwave = (2 * (iwave + off - 1)) * a + b
		flux[iwave] = flux[iwave] / (dwave * area)
	    }

	} else {
	    do iwave = 1, nwave {
		dwave = (2 * (iwave + off - 1)) * a + b
		flux[iwave] = 10.0 ** (-0.4 * flux[iwave]) / (dwave * area)
	    }
	}
end

# PHOTTOCNT -- Special conversion routine for photlam to count units

procedure phottocnt (form, coef, nwave, wave, flux)

char	form[ARB]	# i: output form
char	coef[ARB]	# i: wavelength set coefficients
int	nwave		# i: length of wavelength and flux arrays
real	wave[ARB]	# i: array of wavelengths
real	flux[ARB]	# u: array of fluxes
#--
int	iwave, mwave, off
real	a, b, c, area, dwave

int	stepoffset(), is_magunit()

begin
	# The coefficients define a quadratic function used to
	# define the wavelength set. The wavelength spacing is
	# computed from the derivative of the quadratic

	call gethstarea (area)
	call calcstep (coef, a, b, c, mwave)
	off = stepoffset (wave[1], a, b, c)

	# If the form is not in magnitude units, it is COUNTS
	# If it is, it is OBMAG

	if (is_magunit (form) == NO) {
	    do iwave = 1, nwave {
		dwave = (2 * (iwave + off - 1)) * a + b
		flux[iwave] = flux[iwave] * dwave * area
	    }

	} else {
	    do iwave = 1, nwave {
		dwave = (2 * (iwave + off - 1)) * a + b

		if (flux[iwave] <= 0.0) {
		    flux[iwave] = 100.0
		} else {
		    flux[iwave] = -2.5 * alog10(flux[iwave] * dwave * area)
		}
	    }
	}
end

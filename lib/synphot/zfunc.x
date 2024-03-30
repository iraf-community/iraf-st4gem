# ZFUNC -- Redshift a spectrum by a specified Z value

procedure zfunc (inspec, z, nwave, wave, spec)

real	inspec[ARB]	# i: spectrum to be redshifted
real	z		# i: z value to shift it by
int	nwave		# i: length of wavelength and spectral arrays
real	wave[ARB]	# i: wavelengths at which spectrum is computed
real	spec[ARB]	# o: normalized spectrum, in photlam
#--
int	iwave, jwave
double	shiftwave

real	oneinterp()

begin
	# The code divides through by the z value in order to get 
	# the wavelength whose flux is shifted to the current wavelength

	# Jwave is the index to a wavelength just above the wavelength
	# interpolated at. It is used to simplify search for the next
	# value

	jwave = 1
	do iwave = 1, nwave {
	    shiftwave = wave[iwave] / (1.0 + z)
	    spec[iwave] = oneinterp (jwave, shiftwave, nwave, wave, inspec)
	}

end

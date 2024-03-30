#* HISTORY *
#* B.Simon	04-Oct-94	original

# MONOLAM -- Compute the effective monochromatic flux of a passband

# The effective monochromatic flux is defined by the formula:
#
#   PHOTFEM = PHOTFLAM*RWIDTH*(AVGLAM/lambda)*(TPEAK/T(lambda))

procedure monolam (nwave, wave, filt, lambda, tlambda, emflux)

int     nwave           # i: number of wavelengths
real    wave[ARB]       # i: wavelength array
real    filt[ARB]       # i: throughput array
real	lambda		# u: reference wavelength
real	tlambda		# o: throughput at reference wavelength
real	emflux		# o: effective monochromatic flux
#--
real	area, uflux, peak, ewidth, rwidth, avgwav

real	funit(), peaklam(), avglam()

begin
	# Compute unit stimulus

	call get_hstarea (area)
	uflux = funit (area, nwave, wave, filt)

	# Compute peak wavelength

	peak = peaklam (nwave, filt)

	# Compute rectangular width of passband

	call widthlam (nwave, wave, filt, ewidth, rwidth)

	# Compute average wavelength

	avgwav = avglam (nwave, wave, filt)
	if (IS_INDEFR(lambda))
	    lambda = avgwav

	# Compute throughput at input wavelength

	call syninterp (nwave, wave, filt, 1, lambda, tlambda)

	# Compute effective monochromatic flux

	if (lambda <= 0.0 || tlambda <= 0.0) {
	    emflux = 0.0
	} else {
	    emflux = uflux * rwidth * (avgwav / lambda) * (peak / tlambda)
	}

end

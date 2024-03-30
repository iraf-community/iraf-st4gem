#* HISTORY *
#* B. Simon	22-Sep-94	original

# SPECSCALE -- Compute the scale factor and chi squared for a grid spectrum

procedure specscale (nwave, wave, data, weight, flux, scale, chi2)

int	nwave		# i: length of wave, data, and weight arrays
real	wave[ARB]	# i: wavelength array
real	data[ARB]	# i: input spectrum data
real	weight[ARB]	# i: data point weights
real	flux[ARB]	# i: flux from grid spectrum
real	scale		# o: scaling factor for grid spectrum
real	chi2		# o: chi squared error of fit
#--
pointer	sp, sflux
real	dstim,fstim

real	syntegral(), getchisq()

begin
	# Allocate temporary array

	call smark (sp)
	call salloc (sflux, nwave, TY_REAL)

	# Compute scaling factor

	dstim = syntegral (nwave, wave, data)
	fstim = syntegral (nwave, wave, flux)

	# Compute scaled flux array

	if (fstim <= 0.0) {
	    scale = INDEFR
	    call aclrr (Memr[sflux], nwave)
	} else {
	    scale = dstim / fstim
	    call amulkr (flux, scale, Memr[sflux], nwave)
	}

	# Compute squared difference between data and scaled flux

	chi2 = getchisq (nwave, data, weight, Memr[sflux])

	call sfree (sp)
end


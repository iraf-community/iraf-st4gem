#* HISTORY *
#* B.Simon	12-Jul-94	original

# FWHMLAM2 -- Compute equivalent full width, half max for two bandpasses

real procedure fwhmlam2 (nwave, wave, nband, band1, band2)

int	nwave		# i: length of wavelength set
real	wave[ARB]	# i: wavelength set
int	nband		# i: Number of bandpasses
real	band1[ARB]	# i: First bandpass 
real	band2[ARB]	# i: Second bandpass (or zero)
#--
pointer	sp, band
real	result

real	fwhmlam()

begin
	if (nband == 1) {
	    result = fwhmlam (nwave, wave, band1)

	} else {
	    call smark (sp)
	    call salloc (band, nwave, TY_REAL)

	    call asubr (band1, band2, Memr[band], nwave)
	    result = fwhmlam (nwave, wave, Memr[band])

	    call sfree (sp)
	}

	return (result)
end

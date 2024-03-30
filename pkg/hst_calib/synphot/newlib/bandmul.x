# BANDMUL -- Compute the product of a bandpass and a spectrum

procedure bandmul (nwave, nband, band1, band2, flux)

int	nwave		# i: Length of bandpass and flux arrays
int	nband		# i: Number of bandpasses
real	band1[ARB]	# i: First bandpass
real	band2[ARB]	# i: Second bandpass
real	flux[ARB]	# u: Spectral flux
#--
pointer	sp, flux2

real    tozero()
extern	tozero

begin
	call smark (sp)
	call salloc (flux2, nwave, TY_REAL)

	# If there are two bandpasses, the result is the ratio of the first
	# bandpass times the flux divided by the second bandpass times the
	# flux

	if (nband == 1) {
	    call amulr (band1, flux, flux, nwave)

	} else {
	    call amulr (band2, flux, Memr[flux2], nwave)
	    call amulr (band1, flux, flux, nwave)
	    call advzr (flux, Memr[flux2], flux, nwave, tozero)
	}

	call sfree (sp)
end

# TOZERO - Set error cases to zero

real procedure tozero (val)

real	val		# numerator in division
#--

begin
	return (0.0)
end

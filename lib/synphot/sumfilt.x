# SUMFILT -- Compute numerical passband integral
#
# SUMFILT = Sum [ FILT(I) * WAVE(I) ** NPOW * DWAVE(I) ]

real procedure sumfilt (nwave, wave, npow, band)
 
int	nwave		# i: number of array elements
real	wave[ARB]	# i: input wavelength set array
int	npow		# i: exponential power of wavelength
real	band[ARB]	# i: input passband transmission array
#--
pointer	sp, prod
real	sum

real	syntegral()

begin
	# Allocate array to hold intermediate result

	call smark (sp)
	call salloc (prod, nwave, TY_REAL)

	# Compute appropriate intermediate result 
	# and call integration routine

	switch (npow) {
	case -1:
	    call adivr (band, wave, Memr[prod], nwave)
	    sum = syntegral (nwave, wave, Memr[prod])

	case 0:
	    sum = syntegral (nwave, wave, band)

	case 1:
	    call amulr (band, wave, Memr[prod], nwave)
	    sum = syntegral (nwave, wave, Memr[prod])

	default:
	    call apowkr (wave, npow, Memr[prod], nwave)
	    call amulr (band, Memr[prod], Memr[prod], nwave)
	    sum = syntegral (nwave, wave, Memr[prod])
	}

	call sfree (sp)
	return (sum)

end

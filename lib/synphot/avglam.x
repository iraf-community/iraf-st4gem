# AVGLAM -- Compute the average wavelength of a passband

#  The average wavelength of a bandpass is given by the formula:
#  avglam = Int lam P(lam) dlam / Int P(lam) dlam
 
real procedure avglam (nwave, wave, thruput)

int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavelength array
real	thruput[ARB]	# i: throughput array
#--
pointer	sp, fun
real	sum, sum2, avg

real	syntegral()

begin
	call smark (sp)
	call salloc (fun, nwave, TY_REAL)

	call aabsr (thruput, Memr[fun], nwave)
	sum = syntegral (nwave, wave, Memr[fun])

	call amulr (Memr[fun], wave, Memr[fun], nwave)
	sum2 = syntegral (nwave, wave, Memr[fun])

	if (sum <= 0.0) {
	    avg = 0.0
	} else {
	    avg = sum2 / sum
	}

	call sfree (sp)
	return (avg)
end

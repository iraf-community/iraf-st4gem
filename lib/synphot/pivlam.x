# PIVLAM -- Compute the pivot wavelength of a passband

#  The pivot wavelength should be used to convert fluxes between fnu 
#  and flam. The pivot wavelength is given by the formula
#  the formula:
#
#  pivlam = sqrt( INT(lam * thru * dlam) / INT(thru * dlam / lam) )
 
real procedure pivlam (nwave, wave, thruput)

int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavelength array
real	thruput[ARB]	# i: throughput array
#--
pointer	sp, fun, fun2
real	sum, sum2, pivot

real	syntegral()

begin
	call smark (sp)
	call salloc (fun, nwave, TY_REAL)
	call salloc (fun2, nwave, TY_REAL)

	call aabsr (thruput, Memr[fun], nwave)

	call adivr (Memr[fun], wave, Memr[fun2], nwave)
	sum = syntegral (nwave, wave, Memr[fun2])

	call amulr (Memr[fun], wave, Memr[fun2], nwave)
	sum2 = syntegral (nwave, wave, Memr[fun2])

	if (sum <= 0.0) {
	    pivot = 0.0
	} else {
	    pivot = sqrt (sum2 / sum)
	}

	call sfree (sp)
	return (pivot)
end

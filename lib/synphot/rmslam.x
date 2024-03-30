# RMSLAM -- Compute the rms bandwidth of a passband
#
# The definition used is that from the WF/PC handbook,
# which corresponds to the standard deviation of the passband
# as a function of ln{lam}.
#
# The rms bandwidth is defined in the WF/PC instrument handbook as
#       rmslam^2 = barlam^2 <[ln{lam/barlam}]^2> ,
# where barlam is the bar wavelength, defined by 
# Schneider, Gunn and Hoessel ( 1983 ApJ 264,337 ) as
#       barlam = exp{<ln{lam}>} ,
# and < > denotes an average over the passband in ln{lam}
# <X{lam}> = ( Int P{lam} X{lam} dlam/lam ) / ( Int P{lam} dlam/lam ) .

real procedure rmslam (nwave, wave, thruput)

int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavelength array
real	thruput[ARB]	# i: throughput array
#--
int	iwave
pointer	sp, fun, fun2
real	sumthru, barlam, var

real	syntegral()

begin
	call smark (sp)
	call salloc (fun, nwave, TY_REAL)
	call salloc (fun2, nwave, TY_REAL)

	# Compute total throughput for expectations

	call aabsr (thruput, Memr[fun], nwave)
	call adivr (Memr[fun], wave, Memr[fun], nwave)

	sumthru = syntegral (nwave, wave, Memr[fun])

	if (sumthru <= 0.0) {
	    call sfree (sp)
	    return (0.0)
	}

	# Compute barlam = exp(< ln(lam) >)

	do iwave = 1, nwave 
	    Memr[fun2+iwave-1] = Memr[fun+iwave-1] * log (wave[iwave])

	barlam = syntegral (nwave, wave, Memr[fun2])
	barlam = exp (barlam / sumthru)

	# Compute rms bandwidth ^ 2 = barlam^2 * < ln(lam/barlam)^2 >

	do iwave = 1, nwave 
	    Memr[fun2+iwave-1] = Memr[fun+iwave-1] * 
				log(wave[iwave] / barlam) ** 2 

	var = syntegral (nwave, wave, Memr[fun2]) 
	var = barlam ** 2 * (var / sumthru)

	call sfree (sp)
	return (sqrt(var))
end

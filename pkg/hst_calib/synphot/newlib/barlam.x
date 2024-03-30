#* HISTORY *
#* B.Simon	19-May-94	Adapted from rmslam

# BARLAM -- Computes bar wavelength of a broad bandpass
#
# Barlam is the bar wavelength, defined by Schneider, Gunn and Hoessel 
# (1983 ApJ 264,337) as barlam = exp{<ln{lam}>}, where < > denotes an 
# average over the passband in ln{lam}
# <X{lam}> = ( Int P{lam} X{lam} dlam/lam ) / ( Int P{lam} dlam/lam ) .

real procedure barlam (nwave, wave, filt)

int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavelength array
real	filt[ARB]	# i: throughput array
#--
int	iwave
pointer	sp, fun, fun2
real	sum1, sum2, result

real	syntegral()

begin
	call smark (sp)
	call salloc (fun, nwave, TY_REAL)
	call salloc (fun2, nwave, TY_REAL)

	# Compute total throughput for expectations

	call aabsr (filt, Memr[fun], nwave)
	call adivr (Memr[fun], wave, Memr[fun], nwave)

	sum1 = syntegral (nwave, wave, Memr[fun])

	# Compute barlam = exp(< ln(lam) >)

	do iwave = 1, nwave 
	    Memr[fun2+iwave-1] = Memr[fun+iwave-1] * log (wave[iwave])

	sum2 = syntegral (nwave, wave, Memr[fun2])

	if (sum1 > 0.0) {
	    result = exp (sum2 / sum1)
	} else {
	    result = 0.0
	}

	call sfree (sp)
	return (result)
end

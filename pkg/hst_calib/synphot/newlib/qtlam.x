#* HISTORY *
#* B.Simon	20-Jul-94	original

# QTLAM -- Compute the dimensionless efficiency of a passband
#
# Qtlam is the dimensionless efficiency, as defined in Chapter 6 of the
# WFPC2 Instrument Handbook on p.55. It is defined by the formula
#            qtlam = INT {P{lam} * dlam / lam}

real procedure qtlam (nwave, wave, filt)

int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavelength array
real	filt[ARB]	# i: throughput array
#--
pointer	sp, fun
real	result

real	syntegral()

begin
	call smark (sp)
	call salloc (fun, nwave, TY_REAL)

	# Compute dimensionless efficiency

	call aabsr (filt, Memr[fun], nwave)
	call adivr (Memr[fun], wave, Memr[fun], nwave)

	result = syntegral (nwave, wave, Memr[fun])

	call sfree (sp)
	return (result)
end

#* HISTORY *
#* B.Simon	17-May-94	Original
#* V.Laidler     9-Apr-07	Expanded comments with details & warning

# EFFLAM -- Compute the effective wavelength of a filter and spectrum
#
# The effective wavelength is defined by the formula:
#
# EFFLAM = INT [LAM * FILT * SPEC * LAM * DLAM] /
#          INT [FILT * SPEC * LAM * DLAM]
#
# This formula for the effective wavelength is taken from page
# 836 of the Koornneef et al 1986 paper (Highlights of Astronomy, 1986).
######################################################################
# *** THIS FORMULA IS INCORRECT IF THE FLUX IS IN PHOTLAM ***
# It is correct ONLY if the units of the flux are ergs: that is, if the
# form is "flam".
#
# The calling routine must first convert the flux into the correct
# units!
######################################################################

real procedure efflam (nwave, wave, filt, spec)

int	nwave		# i: number of wavelengths in wave
real	wave[ARB]	# i: wavelength array
real	filt[ARB]	# i: Instrument mode filter array
real	spec[ARB]	# i: spectrum
#--
pointer	sp, prod
real	sum1, sum2, result

real	sumfilt()

begin
	call smark (sp)
	call salloc (prod, nwave, TY_REAL)

	call aabsr (filt, Memr[prod], nwave)
	call amulr (spec, Memr[prod], Memr[prod], nwave)

	sum1 = sumfilt (nwave, wave, 2, Memr[prod])
	sum2 = sumfilt (nwave, wave, 1, Memr[prod])

	if (sum2 > 0.0) {
	    result = sum1 /sum2
	} else {
	    result = 0.0
	}

	call sfree (sp)
	return (result)
end

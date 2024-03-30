#* HISTORY *
#* B.Simon	15-Jun-94	original

# SETPHOTFORM -- Convert effective stimulus to the requested output form

procedure setphotform (inform, outform, wave, nwave, band1, band2, nband, 
		       stim, nstim)

char	inform[ARB]	# i: input form
char	outform[ARB]	# i: output form
real	wave[ARB]	# i: wavelength array
int	nwave		# i: number of wavelengths
real	band1[ARB]	# i: first bandpass
real	band2[ARB]	# i: second bandpass (or zero)
int	nband		# i: number of bandpasses
real	stim[ARB]	# u: effective stimulus
int	nstim		# i: length of stimulus array
#--
bool	noset
int	istim
pointer	sp, form1, form2
real	irate, orate, scale, term

bool	streq()
int	is_magunit()
real	photrate()

begin
	# If the input and output forms are the same,
	# we can skip the conversion

	call smark (sp)
	call salloc (form1, SZ_FNAME, TY_CHAR)
	call salloc (form2, SZ_FNAME, TY_CHAR)

	call strcpy (inform, Memc[form1], SZ_FNAME)
	call strcpy (outform, Memc[form2], SZ_FNAME)

	call strfix (Memc[form1])
	call strfix (Memc[form2])

	noset = streq (Memc[form1], Memc[form2])
	call sfree (sp)

	if (noset)
	    return

	# Compute the photometric rate for a unit flux
	# The ratio of these rates is the scaling factor 

	irate = photrate (inform, wave, nwave, band1, band2, nband)
	orate = photrate (outform, wave, nwave, band1, band2, nband)

	# Compute the scaling factor and apply to the stimulus 
	# to convert its form

	if (orate == 0.0) {
	    call amovkr (INDEFR, stim, nstim)

	} else if (irate == 0.0) {
	    if (is_magunit (outform) == NO) {
		call amovkr (0.0, stim, nstim)
	    } else {
		call amovkr (100.0, stim, nstim)
	    }

	} else {
	    scale = irate / orate
	    term = -2.5 * alog10 (scale)

	    if (is_magunit (inform) == NO) {
		if (is_magunit (outform) == NO) {
		    do istim = 1, nstim
			stim[istim] = stim[istim] * scale
		} else {
		    do istim = 1, nstim {
			if (stim[istim] > 0.0) {
			    stim[istim] = -2.5 * alog10 (stim[istim]) + term
			} else {
			    stim[istim] = 100.0
			}
		    }
		}

	    } else {
		if (is_magunit (outform) == NO) {
		    do istim = 1, nstim
			stim[istim] = 10.0 ** (-0.4 * (stim[istim] + term))
		} else {
		    do istim = 1, nstim
			stim[istim] = stim[istim] + term
		}
	    }
	}
end

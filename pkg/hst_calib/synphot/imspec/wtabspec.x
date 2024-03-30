#* HISTORY *
#* B.Simon	03-Mar-93	original
#* B.Simon	07-Jan-94	add output form
#* B.Simon	03-Jun-94	revised for new synphot library
#* B.Simon	21-Aug-96	calls rebin

# WTABSPEC -- Write the spectrum to an sdas table

procedure wtabspec (table, outform, olength, ilen, iwave, ispec)

char	table[ARB]	# i: table name
char	outform[ARB]	# i: output spectral form string
int	olength		# i: length of output spectrum
int	ilen		# i: number of wavelengths
pointer	iwave		# i: wavelength array
pointer	ispec		# i: spectral array
#--
int	opix, olen, done
pointer	sp, owave, ospec, tp, wavptr, fluxptr
real	crval, cd

string	blank    " "

int	phottoany()
pointer	tbtopn()

begin
	call smark (sp)

	# Compute length of output file

	if (IS_INDEFI(olength)) {
	    olen = ilen
	} else {
	    olen = olength
	}

	# Compute new wavelength scale if output length different
	# than input length. Then interpolate on new scale.

	if (ilen == olen) {
	    owave = iwave
	    ospec = ispec

	} else {
	    call salloc (owave, olen, TY_REAL)
	    call salloc (ospec, olen, TY_REAL)

	    crval = Memr[iwave]
	    cd = (Memr[iwave+ilen-1] - Memr[iwave]) / real (olen - 1)

	    do opix = 0, olen-1
		Memr[owave+opix] = crval + cd * real (opix)

	    call rebin (ilen, Memr[iwave], Memr[ispec], 
			olen, Memr[owave], Memr[ospec])
	}

	# Convert spectrum to output units

	done = phottoany (outform, olen, Memr[owave], Memr[ospec])

	# Create output table and its columns

	tp = tbtopn (table, NEW_FILE, 0)
	call tbcdef (tp, wavptr, "WAVELENGTH", "angstroms", blank, 
		     TY_REAL, 1, 1)
	call tbcdef (tp, fluxptr, "FLUX", outform, blank, TY_REAL, 1, 1)
	call tbtcre (tp)

	# Write spectrum to table

	call tbcptr (tp, wavptr, Memr[owave], 1, olen)
	call tbcptr (tp, fluxptr, Memr[ospec], 1, olen)

	# Close table and free memory

	call tbtclo (tp)
	call sfree (sp)
end

#* HISTORY *
#* B.Simon	16-May-94	Original
#* B.Simon	20-Jul-94	added qtlam
#* V.Laidler     9-Apr-07       Replace EFFLAM by EFFLERG, EFFLPHOT
#*
# PHOTSTIM -- Compute the effective stimulus for calcphot

procedure photstim (form, func, nwave, wave, nband, band1, band2, flux, result)
		    
char	form[ARB]	# i: Form of result
char	func[ARB]	# i: Output function
int	nwave		# i: length of wavelength set
real	wave[ARB]	# i: wavelength set
int	nband		# i: Number of bandpasses
real	band1[ARB]	# i: First bandpass 
real	band2[ARB]	# i: Second bandpass (or zero)
real	flux		# i: Spectral flux
real	result		# o: Computed result
#--
int	done
pointer	sp, band, prod

string	badfunc  "photstim: illegal function"
string    flamform "flam"

int	phottoany(), findfunc()
real	rmslam(), fwhmlam(), barlam(), avglam()
real	efflam(), effstim2()

begin
	call smark (sp)
	call salloc (band, nwave, TY_REAL)
	call salloc (prod, nwave, TY_REAL)

	# Multiply bandpass by the flux

	if (nband == 1) {
	    call amovr (band1, Memr[band], nwave)
	} else {
	    call asubr (band1, band2, Memr[band], nwave)
	}

	call amovr (flux, Memr[prod], nwave)
	done = phottoany (form, nwave, wave, Memr[prod])
	call amulr (Memr[band], Memr[prod], Memr[prod], nwave)

	# Compute user selected output

	switch (findfunc (func)) {
	case 1: # effstim
	    result = effstim2 (form, nwave, wave, nband, 
			       band1, band2, flux)
	case 2: # rmslam
	    result = rmslam (nwave, wave, Memr[prod])
	case 3: # fwhmlam
	    result = fwhmlam (nwave, wave, Memr[prod])
	case 4: # barlam
	    result = barlam (nwave, wave, Memr[prod])
	case 5: # avglam
	    result = avglam (nwave, wave, Memr[prod])
	case 6: # old efflam = calculated with flux in photlams
	    result = efflam (nwave, wave, Memr[band], flux)
	case 7: # efflam as per Koornneef et al 1986 p 836
	    call amovr (flux, Memr[prod], nwave)
	    done = phottoany (flamform, nwave, wave, Memr[prod])
	    result = efflam (nwave, wave, Memr[band], Memr[prod])
	default: # oops!!
	    call printerr_str (badfunc, func)
	}

	call sfree (sp)
end


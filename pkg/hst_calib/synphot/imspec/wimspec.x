#* HISTORY *
#* B.Simon	05-Mar-93	original
#* B.Simon	07-Jan-94	add output form
#* B.Simon	03-Jun-94	revised for new synphot library
#* B.Simon	26-Oct-94	changed CTYPE value written to image
#* B.Simon	21-Aug-96	calls rebin
#* H.Bushouse	15-May-01	changed CTYPE value from "LAMBDA" to "LINEAR"

# WIMSPEC -- Write the spectrum to an image

procedure wimspec (image, wave, outform, olength, ilen, iwave, ispec)
		   

char	image[ARB]	# i: image name
char	wave[ARB]	# i: optional wavelength image
char	outform[ARB]	# i: spectral form string
int	olength		# i: length of output spectrum
int	ilen		# i: number of wavelengths
pointer	iwave		# i: wavelength array
pointer	ispec		# i: spectral array
#--
bool	linear
int	opix, olen, done
pointer	sp, owave, ospec
real	crval, cd

bool	is_linear()
int	phottoany()

begin
	# Allocate dynamic memory

	call smark (sp)

	# Compute length of output file

	if (IS_INDEFI(olength)) {
	    olen = ilen
	} else {
	    olen = olength
	}

	# Compute new wavelength scale if output length different
	# than input length. Then interpolate on new scale.

	crval = Memr[iwave]
	cd = (Memr[iwave+ilen-1] - Memr[iwave]) / real (olen - 1)

	linear = is_linear (Memr[iwave], ilen)

	if (wave[1] != EOS || (linear && ilen == olen)) {
	    owave = iwave
	    ospec = ispec

	} else {
	    call salloc (owave, olen, TY_REAL)
	    call salloc (ospec, olen, TY_REAL)

	    do opix = 0, olen-1
		Memr[owave+opix] = crval + cd * double(opix)

	    call rebin (ilen, Memr[iwave], Memr[ispec], 
			olen, Memr[owave], Memr[ospec])
	}

	# Convert spectrum to output units

	done = phottoany (outform, olen, Memr[owave], Memr[ospec])

	# Write either one or two images depending on whether
	# wavelength file is present

	if (wave[1] == EOS) {
	    call wrtonedim (image, crval, cd, "LINEAR", Memr[ospec], olen)

	} else {
	    call wrtonedim (wave, 1.0, 1.0, "PIXEL", Memr[owave], olen)
	    call wrtonedim (image, 1.0, 1.0, "PIXEL", Memr[ospec], olen)
	}

	call sfree (sp)
end

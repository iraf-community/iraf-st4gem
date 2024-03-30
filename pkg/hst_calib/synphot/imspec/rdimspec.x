include	<imhdr.h>

#* HISTORY*
#* B.Simon	04-Mar-93	original
#* B.Simon	07-Jan-94	add input form
#* B.Simon	03-Jun-94	revised for new synphot library
#* B.Simon	26-Oct-94	removed wavelength units conversion

# RDIMSPEC -- Read a spectrum from an image

procedure rdimspec (image, wave, inform, badpix, ilen, iwave, ispec)

char	image[ARB]	# i: image name
char	wave[ARB]	# i: optional wavelength image
char	inform[ARB]	# i: spectral form string
real	badpix		# i: replacement value for bad pixel
int	ilen		# o: number of wavelengths
pointer	iwave		# o: wavelength array
pointer	ispec		# o: spectral array
#--
int	ipix, done
double	crpix, crval, cd
pointer	sp, ctype, imi, imw, ibuf, wbuf

string	baddimen  "Image is not one-dimensional"
string	badlength "Input and wavelength images are not same length"

int	checkdim(), fillnull(), anytophot()
pointer	immap(), imgl1r()

begin
	call smark (sp)
	call salloc (ctype, SZ_FNAME, TY_CHAR)

	# Read the image containing the flux

	imi = immap (image, READ_ONLY, 0)
	if (checkdim (imi) > 1) {
	    call imunmap (imi)
	    call printerr_str (baddimen, image)
	}

	ilen = IM_LEN(imi, 1)
	call malloc (iwave, ilen, TY_REAL)
	call malloc (ispec, ilen, TY_REAL)

	ibuf = imgl1r (imi)
	call amovr (Memr[ibuf], Memr[ispec], ilen)

	# If no wavelength file, calculate wavelengths from cd matrix
	# Otherwise read them from the wavelength file

	if (wave[1] == EOS) {
	    call stx_getcoord (imi, crpix, crval, cd, 1, Memc[ctype], SZ_FNAME)

	    do ipix = 1, ilen
		Memr[iwave+ipix-1] = crval + (double(ipix) - crpix) * cd

	} else {
	    imw = immap (wave, READ_ONLY, 0)

	    if (checkdim (imw) > 1) {
		call imunmap (imi)
		call imunmap (imw)
		call mfree (iwave, TY_REAL)
		call mfree (ispec, TY_REAL)

		call printerr_str (baddimen, wave)
	    }

	    if (ilen != IM_LEN(imw,1)) {
		call imunmap (imi)
		call imunmap (imw)
		call mfree (iwave, TY_REAL)
		call mfree (ispec, TY_REAL)

		call printerr_str (badlength, wave)
	    }

	    Memc[ctype] = EOS
	    wbuf = imgl1r (imw)

	    call amovr (Memr[wbuf], Memr[iwave], ilen)
	    call imunmap (imw)
	}

	# Convert the units of the image to photlam

	done = fillnull (badpix, ilen, Memr[iwave], Memr[ispec])
	done = anytophot (inform, ilen, Memr[iwave], Memr[ispec])

	# Put the spectrum in ascending order of wavelength

	call synsort2 (ilen, Memr[iwave], Memr[ispec])

	call imunmap (imi)
	call sfree (sp)

end

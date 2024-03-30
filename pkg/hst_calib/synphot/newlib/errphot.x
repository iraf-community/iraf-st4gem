#* HISTORY *
#* B.Simon	13-Jun-94	original
#* B.Simon	08-Sep-94	take absolute value of error for mag units

# ERRPHOT -- Convert flux and its error to photlam units

procedure errphot (form, nwave, wave, flux, error)

char	form		# i: input flux units
int	nwave		# i: length of wavelength and other arrays
real	wave[ARB]	# i: wavelength array
real	flux[ARB]	# u: flux array
real	error[ARB]	# u: flux error array
#--
int	done
pointer	sp, fluxp, fluxm

int	anytophot()

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (fluxp, nwave, TY_REAL)
	call salloc (fluxm, nwave, TY_REAL)

	# A missing error array is signalled by an INDEF value

	if (IS_INDEFR (error[1])) {
	    done = anytophot (form, nwave, wave, flux)

	} else {
	    # We approximate the conversion of error to photlam units
	    # by adding and subtracting the error from the flux, converting
	    # these arrays to photlam, and setting the converted error to
	    # half the difference of the converted arrays. This is an 
	    # approximation that is valid when the error is small relative
	    # to the flux.

	    call aaddr (flux, error, Memr[fluxp], nwave)
	    call asubr (flux, error, Memr[fluxm], nwave)

	    done = anytophot (form, nwave, wave, flux)
	    done = anytophot (form, nwave, wave, Memr[fluxp])
	    done = anytophot (form, nwave, wave, Memr[fluxm])

	    call asubr (Memr[fluxp], Memr[fluxm], error, nwave)
	    call amulkr (error, 0.5, error, nwave)
	    call aabsr (error, error, nwave)
	}

	call sfree (sp)
end

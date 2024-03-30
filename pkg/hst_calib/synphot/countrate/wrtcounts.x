#* HISTORY *
#* B.Simon	31-May-94	original

# WRT_COUNTS -- Write results from countrate

procedure wrt_counts (mode, spectrum, form, refwave, 
		      pivot, fwhm, flux_tot, flux_ref)

char	mode[ARB]	# i: observation mode
char	spectrum[ARB]	# i: spectrum
char	form[ARB]	# i: output form
real	refwave		# i: reference wavelength
real	pivot		# i: pivot wavelength
real	fwhm		# i: full width half maximum
real	flux_tot	# i: total flux
real	flux_ref	# i: flux at reference wavelength
#--
pointer	sp, uform

string	hdrfmt1   "   Pivot       Equiv Gaussian  Total Flux\n"
string	hdrfmt2   " Wavelength         FWHM         %s\n"
string	hdrfmt3   "   Pivot       Equiv Gaussian  Total Flux  Flux at %d A\n"
string	hdrfmt4   " Wavelength         FWHM         %s        %s\n"
string	datafmt1  "%10g    %12g     %10g\n"
string	datafmt2  "%10g    %12g     %10g    %10g\n"

begin
	# Convert form to upper case for output

	call smark (sp)
	call salloc (uform, SZ_FNAME, TY_CHAR)

	call strcpy (form, Memc[uform], SZ_FNAME)
	call strupr (Memc[uform])

	# Print observation mode and spectrum

	call printf ("Mode = %s\n")
	call pargstr (mode)

	call printf ("Spectrum:  %s\n")
	call pargstr (spectrum)

	# Print results

	if (IS_INDEFR(refwave)) {
	    call printf (hdrfmt1)
	    call printf (hdrfmt2)
	    call pargstr (Memc[uform])

	    call printf (datafmt1)
	    call pargr (pivot)
	    call pargr (fwhm)
	    call pargr (flux_tot)

	} else {
	    call printf (hdrfmt3)
	    call pargr (refwave)

	    call printf (hdrfmt4)
	    call pargstr (Memc[uform])
	    call pargstr (Memc[uform])

	    call printf (datafmt2)
	    call pargr (pivot)
	    call pargr (fwhm)
	    call pargr (flux_tot)
	    call pargr (flux_ref)
	}

	call sfree (sp)
end

#* HISTORY *
#* B.Simon	20-Sep-00	Updated calculation of wavelength set

# PHOTBAND -- Calculate and print the bandpass input to calcphot

procedure photband (grtbl, cmptbl, mode, nwave, wave, maxcode, nband, 
		    mode1, mode2, ocode1, ocode2, band1, band2)

char	grtbl[ARB]	# i: Instrument graph table
char	cmptbl[ARB]	# i: Component lookup table
char	mode[ARB]	# i: Observation mode expression
int	nwave		# i: length of wavelength set
real	wave[ARB]	# i: wavelength set
int	maxcode		# i: Declared length of code arrays
int	nband		# o: Number of passbands
char	mode1[ARB]	# o: First part of obsmode expression
char	mode2[ARB]	# o: Second part of obsmode expression
int	ocode1[ARB]	# o: Code for first obsmode expression
int	ocode2[ARB]	# o: Code for second obsmode expression (or NULL)
real	band1[ARB]	# o: First bandpass 
real	band2[ARB]	# o: Second bandpass (or zero)
#--
int	degree
real	pivot1, fwhm1, pivot2, fwhm2

string	notband  "Not a passband"
string	mode1fmt "Mode = %s\n"
string	mode2fmt "Mode = %s - %s\n"
string	hdrfmt	 "   Pivot       Equiv Gaussian\n Wavelength         FWHM\n"
string	datafmt  "%10g    %12g    %s\n"

real	pivlam(), fwhmlam()

begin
	# Split obsmode at subtraction sign so that
	# two halves can be computed separately

	call splitexp (mode, nband, mode1, mode2, maxcode)
			   
	call syncompile (mode1, ocode1, maxcode)
	if (nband == 1) {
	    ocode2[1] = -1
	} else {
	    call syncompile (mode2, ocode2, maxcode)
	}

	# Calculate passband, its pivot wavelength, and fwhm
	# for the first passband

	call syncalc (ocode1, maxcode, NULL, nwave, wave, grtbl, cmptbl,
		      band1, degree)

	if (degree != 0)
	    call printerr_str (notband, mode1)

	pivot1 = pivlam (nwave, wave, band1)
	fwhm1 = fwhmlam (nwave, wave, band1)

	if (nband == 1) {
	    # Print results if only one observation mode

	    call amovkr (0.0, band2, nwave)

	    call printf (mode1fmt)
	    call pargstr (mode1)

	    call printf (hdrfmt)

	    call printf (datafmt)
	    call pargr (pivot1)
	    call pargr (fwhm1)
	    call pargstr (mode1)

	} else {
	    # Calculate passband, its pivot wavelength, and fwhm
	    # for the second passband

	    call syncalc (ocode2, maxcode, NULL, nwave, wave, grtbl, cmptbl,
			  band2, degree)

	    if (degree != 0)
		call printerr_str (notband, mode2)

	    pivot2 = pivlam (nwave, wave, band2)
	    fwhm2 = fwhmlam (nwave, wave, band2)

	    # Print results for both observation modes

	    call printf (mode2fmt)
	    call pargstr (mode1)
	    call pargstr (mode2)

	    call printf (hdrfmt)

	    call printf (datafmt)
	    call pargr (pivot1)
	    call pargr (fwhm1)
	    call pargstr (mode1)

	    call printf (datafmt)
	    call pargr (pivot2)
	    call pargr (fwhm2)
	    call pargstr (mode2)
	}

	call flush (STDOUT)
end

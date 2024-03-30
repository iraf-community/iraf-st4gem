#* HISTORY *
#* B.Simon	08-Jun-94	adapted from photband.x

# SYNCALC2 -- Calculate two bandpass expressions

procedure syncalc2 (graftab, comptab, mode1, mode2, ocode1, ocode2, maxcode, 
		    nwave, wave, nband, band1, band2)

char	graftab[ARB]	# i: Instrument graph table
char	comptab[ARB]	# i: Component lookup table
char	mode1[ARB]	# i: First obsmode expression
char	mode2[ARB]	# i: Second obsmode expression
int	ocode1[ARB]	# i: Code for first obsmode expression
int	ocode2[ARB]	# i: Code for second obsmode expression (or NULL)
int	maxcode		# i: Declared length of code arrays
int	nwave		# i: length of wavelength set
real	wave[ARB]	# i: wavelength set
int	nband		# i: Number of obsmode passbands
real	band1[ARB]	# o: First bandpass 
real	band2[ARB]	# o: Second bandpass (or zero)
#--
int	degree

string	notband  "Not a passband"

begin
	# Calculate first passband

	call syncalc (ocode1, maxcode, NULL, nwave, wave, graftab, comptab,
		      band1, degree)
	if (degree != 0)
	    call printerr_str (notband, mode1)

	# Optionally calculate second passband

	if (nband == 1) {
	    call amovkr (0.0, band2, nwave)

	} else{
	    call syncalc (ocode2, maxcode, NULL, nwave, wave, graftab, comptab,
			  band2, degree)
	    if (degree != 0)
		call printerr_str (notband, mode2)
	}

end

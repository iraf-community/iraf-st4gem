# WAVEMAX -- Maximum of two wavelengths, which may be INDEF

real procedure wavemax (w1, w2)

real	w1		# i: First wavelength
real	w2		# i: Second wavelength
#--
real	wmax

begin
	if (IS_INDEFR (w1)) {
	    wmax = w2
	} else if (IS_INDEFR (w2)) {
	    wmax = w1
	} else {
	    wmax = max (w1, w2)
	}

	return (wmax)
end

# WAVEMIN -- Minimum of two wavelengths, which may be INDEF

real procedure wavemin (w1, w2)

real	w1		# i: First wavelength
real	w2		# i: Second wavelength
#--
real	wmin

begin
	if (IS_INDEFR (w1)) {
	    wmin = w2
	} else if (IS_INDEFR (w2)) {
	    wmin = w1
	} else {
	    wmin = min (w1, w2)
	}

	return (wmin)
end

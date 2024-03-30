# WAVESET - Generate an array of wavelengths between user supplied limits

procedure waveset (logspace, minwave, maxwave, nwave, wave)

bool	logspace	# i: use logaraithmic spacing between wavelengths?
real	minwave		# i: minimum wavelenth
real	maxwave		# i: maximum wavelength
int	nwave		# i: number of wavelengths
real	wave[ARB]	# o: wavelength set
#--
int	iwave
real	w1, w2, frac

string	norange  "Computed wavelength range is zero"

errchk	synphoterr

begin

	if (minwave >= maxwave)
	    call synphoterr (norange, "waveset")

	if (nwave == 1) {
	    # Check for pathalogical case
	    wave[1] = 0.5 * (minwave + maxwave)

	} else if (logspace) {
	    # Space wavelengths logarithmically (more at shorter wavelengths)
	    w1 = log10 (minwave)
	    w2 = log10 (maxwave)

	    wave[1] = minwave
	    wave[nwave] = maxwave

	    do iwave = 1, nwave-2 {
		frac = real(iwave) / real(nwave-1)
		wave[iwave+1] = 10 ** (w1 * (1.0 - frac) + w2 * frac)
	    }

	} else {
	    # Space wavelengths linearly (equal distance between wavelengths)
	    wave[1] = minwave
	    wave[nwave] = maxwave

	    do iwave = 1, nwave-2 {
		frac = real(iwave) / real(nwave-1)
		wave[iwave+1] = minwave * (1.0 - frac) + maxwave * frac
	    }
	}

end

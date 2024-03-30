#* HISTORY *
#* B.Simon	06-Aug-96	original

# CROSSWAVE -- Intersect wavelength set with wavelength range

procedure crosswave (minwave, maxwave, nwave, wave)

real	minwave		# i: minimum of wavelength range
real	maxwave		# i: maximum of wavelength range
int	nwave		# u: number of wavelengths
real	wave[ARB]	# u: wavelength set
#--
int	iwave, imin, imax

begin
	# Compute intersection  points

	imin = nwave + 1
	do iwave = 1, nwave {
	    if (wave[iwave] >= minwave) {
		imin = iwave
		break
	    }
	}

	imax = 0
	do iwave = nwave, 1, -1 {
	    if (wave[iwave] <= maxwave) {
		imax = iwave
		break
	    }
	}

	nwave = max (0, (imax - imin) + 1)

	# Copy wavelength set into start of wave array

	if (imin > 1) {
	    do iwave = 1, nwave
		wave[iwave] = wave[imin+iwave-1]
	}

end

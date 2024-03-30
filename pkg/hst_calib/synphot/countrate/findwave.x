#* HISTORY *
#* B.Simon	25-Mar-93	original

# FINDWAVE -- Find the nearest wavelength to the specified wavelength

int procedure findwave (cenwave, wave, nwave)

real	cenwave		# i: wavelength to search for
real	wave[ARB]	# i: sorted array of wavelengths
int	nwave		# i: length of wavelength array
#--
int	iwave
real	oldist, dist

begin
	oldist = abs (cenwave - wave[1])

	for (iwave = 2; iwave <= nwave; iwave = iwave + 1) {
	    dist = abs (cenwave - wave[iwave])
	    if (dist > oldist)
		break

	    oldist = dist
	}

	return (iwave - 1)
end

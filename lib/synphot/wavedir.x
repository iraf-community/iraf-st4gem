define	FAST		NO

# WAVEDIR - Compute direction of wavelength array

int procedure wavedir (nwave, wave)

int	nwave		# i: number of wavelengths
real	wave[ARB]	# i: wavelength array
#--
int	iwave, dir
real	dwave

string	samewave "Duplicate wavelength"
string	sortwave "Wavelength out of sort"

errchk	synrealerr

begin
	# Fast option doesn't check array for monotonicity

	if (FAST == YES) {
	    if (nwave < 2) {
		dir = 0
	    } else {
		dir = sign (1.0, (wave[2] - wave[1]))
	    }

	    return (dir)
	}

	# Check each point in wavelength array to ensure monotonicity

	dir = 0
	do iwave = 2, nwave {
	    dwave = wave[iwave] - wave[iwave-1]
	    switch (dir) {
	    case -1:
		if (dwave == 0.0) {
		    dir = 0
		    call synrealerr (samewave, wave[iwave])
		} else if (dwave > 0.0) {
		    dir = 0
		    call synrealerr (sortwave, wave[iwave])
		}

	    case 0:
		if (dwave == 0.0) {
		    dir = 0
		    call synrealerr (samewave, wave[iwave])
		} else if (dwave > 0.0) {
		    dir = 1
		} else {
		    dir = -1
		}
	    case 1:
		if (dwave == 0.0) {
		    dir = 0
		    call synrealerr (samewave, wave[iwave])
		} else if (dwave < 0.0) {
		    dir = 0
		    call synrealerr (sortwave, wave[iwave])
		}
	    }
	}

	return (dir)
end

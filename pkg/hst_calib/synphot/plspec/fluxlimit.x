define	MAG_FACTOR	15
define	NON_MAG_FACTOR	0.001

#* HISTORY *
#* B.Simon	26-May-00	original
#* B.Simon	07-Jun-00	fixed call to setlimit

# FLUXLIMIT -- Compute plot limits from flux values

procedure fluxlimit (nwave, wave, flux, form, inlimit, outlimit)

int	nwave		# i: Length of wave and spectrum arrays
real	wave[ARB]	# i: Wavelength array
real	flux[ARB]	# i: Flux array
char	form[ARB]	# i: Output flux units
real	inlimit[4]	# i: Input plot limits
real	outlimit[4]	# u: Output plot limis
#--
int	iwave, jwave, kwave, mwave
real	minval, maxval

bool	is_magunit()

begin
	# Find minimum and maximum flux values

	minval = flux[1]
	maxval = flux[1]

	do iwave = 2, nwave {
	    if (flux[iwave] < minval)
		minval = flux[iwave]

	    if (flux[iwave] > maxval)
		maxval = flux[iwave]
	}

	# Set minimum value as fraction of maximum

	if (is_magunit (form)) {
	    if (MAG_FACTOR + minval >= maxval) {
		# No need to truncate, calculate full limits

		jwave = 1
		kwave = nwave

	    } else {
		# Determine where the flux crosses the maximum value

		maxval = MAG_FACTOR + minval
		do iwave = 2, nwave {
		    if (flux[iwave] < maxval) {
			jwave = iwave - 1
			break
		    }
		}

		do iwave = nwave-1, 1, -1 {
		    if (flux[iwave] < maxval) {
			kwave = iwave + 1
			break
		    }
		}
	    }

	} else {
	    if (NON_MAG_FACTOR * maxval <= minval) {
		# No need to truncate, calculate full limits

		jwave = 1
		kwave = nwave

	    } else {
		# Determine where the flux crosses the minimum value

		minval = NON_MAG_FACTOR * maxval
		do iwave = 2, nwave {
		    if (flux[iwave] > minval) {
			jwave = iwave - 1
			break
		    }
		}

		do iwave = nwave-1, 1, -1 {
		    if (flux[iwave] > minval) {
			kwave = iwave + 1
			break
		    }
		}
	    }
	}

	# Call setlimit with a truncated spectrum to calculate
	# limits based on flux values

	mwave = (kwave - jwave) + 1
	call setlimit (mwave, wave[jwave], flux[jwave], inlimit, outlimit)

end

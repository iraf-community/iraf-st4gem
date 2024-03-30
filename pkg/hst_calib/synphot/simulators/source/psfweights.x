#* HISTORY *
#* B.Simon	17-Feb-95	original

# PSFWEIGHTS -- Compute weight factors for interpolating psfs

procedure psfweights (nwave, wave, thruput, flux, numpsf, wavepsf, weight)

int	nwave		# i: length of wavelength and flux arrays
real	wave[ARB]	# i: wavelengths where thruput and flux are evaluated
real	thruput[ARB]	# i: instrument thruput
real	flux[ARB]	# i: object spectral flux
int	numpsf		# i: number of weights
real	wavepsf[ARB]	# i: wavelengths at which psf weights are computed
real	weight[ARB]	# o: weights, scaled to one
#--
int	ipsf, iwave, jwave, nlim
pointer	sp, limit
real	total

string	badweight  "psfweights: bad value for numpsf"

real	effstim()

begin
	# Special case of a single psf sets the weight to one

	if (numpsf <= 0)
	    call printerr_int (badweight, numpsf)

	if (numpsf == 1) {
	    weight[1] = 1.0
	    return
	}

	# Temporary array to hold integration limits

	call smark (sp)
	call salloc (limit, numpsf+1, TY_REAL)

	# Compute integration limits

	Memr[limit] = min (wave[1], wavepsf[1])
	Memr[limit+numpsf] = max (wave[nwave], wavepsf[numpsf])

	do ipsf = 1, numpsf-1
	    Memr[limit+ipsf] = 0.5 * (wavepsf[ipsf+1] + wavepsf[ipsf])

	# The weight of the psf is the flux between integration limits

	ipsf = 1
	iwave = 1
	jwave = 1
	total = 0.0

	while (iwave <= nwave && ipsf <= numpsf) {
	    if (wave[iwave] < Memr[limit+ipsf]) {
		iwave = iwave + 1

	    } else {
		nlim = iwave - jwave
		if (nlim == 0) {
		    weight[ipsf] = 0.0

		} else {
		    weight[ipsf] = effstim (nlim, wave[jwave], thruput[jwave],
					    flux[jwave], "counts")
					    
		    total = total + weight[ipsf]
		    jwave = iwave
		}

		ipsf = ipsf + 1
	    }
	}

	# Set the weight of the psf closest to the wavelength range 
	# to one if all the weights are zero

	if (total == 0.0) {
	    ipsf = min (ipsf, numpsf)
	    weight[ipsf] = 1.0
	    ipsf = ipsf + 1
	    total = 1.0
	}

	# Set remaining weights to zero

	while (ipsf <= numpsf) {
	    weight[ipsf] = 0.0
	    ipsf = ipsf + 1
	}

	# Normalize weights to one

	total = 1.0 / total
	call amulkr (weight, total, weight, numpsf)

	call sfree (sp)
end

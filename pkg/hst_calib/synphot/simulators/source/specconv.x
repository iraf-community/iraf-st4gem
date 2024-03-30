include "grating.h"

#* HISTORY *
#* B.Simon	17-Jul-95	original
#* B.Simon	29-Nov-95	reworked to use dispersion equation
#* B.Simon	24-Jan-96	fixed typos in convolution code

# SPECCONV -- Convolve the spectrum with a line from the source

int procedure specconv (nwave, wave, spec, order, grating, 
			jdx, ndx, ndy, anobj, line)

int	nwave		# i: length of spectrum
real	wave[ARB]	# i: source wavelength set
real	spec[ARB]	# i: source spectrum
int	order		# i: spectral order
real	grating[LEN_G]	# i: grating dispersion coefficients
int	jdx		# i: line in source to convolve
int	ndx		# i: x dimension of source
int	ndy		# i: y dimension of source
real	anobj[ndx,ndy]	# i: masked source 
real	line[ARB]	# o: convolved line
#--
int	iwave, jwave, idy, mdy
pointer	sp, iprof, oprof
real	lineflux, mag

real	magdisp()

begin
	# Calculate line flux
	# Bypass convolution if source line is all zero

	lineflux = 0.0
	do idy = 1, ndy
	    lineflux = lineflux + anobj[jdx,idy]

	if (lineflux <= 0.0)
	    return (NO)

	# The convolution of the spectrum with the object profile
	# in the dispersion direction cannot be done in the 
	# straightforward way. First, we must use the dispersion 
	# equation to compute a magnification factor for the profile.
	# After calculating the magnification factor, we shrink or
	# expand the object profile and convolve with the magnified
	# profile.

	# Compute the magnification factor caused by the dispersion
	# of the grating

	iwave = (nwave + 1)/ 2
	mag = magdisp (ndy, nwave, wave[iwave], order, grating)

	# If the magnified  profile has a length less than one, there
	# is no point in convolving. Simply multiply by the total
	# flux and exit

	mdy = ndy * mag
	if (mdy <= 1) {
	    call amulkr (spec, lineflux, line, nwave)
	    return (YES)
	}

	# Allocate memory

	mdy = mdy + 2

	call smark (sp)
	call salloc (iprof, ndy, TY_REAL)
	call salloc (oprof, mdy, TY_REAL)

	# Copy object profile to one dimensional array

	do idy = 1, ndy {
	    Memr[iprof+idy-1] = anobj[jdx,idy]
	}

	# Resample the profile according to the magnification factor

	call vecrebin (ndy, Memr[iprof], mdy, Memr[oprof])

	# Convolve the spectrum with the resampled profile
	# (One dimensional convolution)

	call aclrr (line, nwave)

	do iwave = 1, nwave {
	    jwave = mdy + (iwave - 1)
	    do idy = 1, mdy {
		if (jwave <= nwave)
		    line[jwave] = line[jwave] + Memr[oprof+idy-1] * 
				  spec[iwave]
		jwave = jwave - 1
	    }
	}

	call sfree (sp)
	return (YES)
end

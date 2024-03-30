#* HISTORY *
#* B.Simon	17-Jul-95	original

# CONVSPEC -- Convolve the spectrum with a line from the source

int procedure convspec (nwave, spec, itx, ntx, nty, asrc, line)

int	nwave		# i: length of spectrum
real	spec[ARB]	# i: source spectrum
int	itx		# i: line in source to convolve
int	ntx		# i: x dimension of source
int	nty		# i: y dimension of source
real	asrc[ntx,nty]	# i: masked source 
real	line[ARB]	# o: convolved line
#--
int	ity, iwave, jwave
real	lineflux, total

real	asumr()

begin
	# Set output array to zero

	call aclrr (line, nwave)

	# Calculate line flux
	# Bypass convolution if source line is all zero

	lineflux = 0.0
	do ity = 1, nty
	    lineflux = lineflux + asrc[itx,ity]

	if (lineflux <= 0.0)
	    return (NO)

	# One dimensional convolution

	do iwave = 1, nwave {
	    jwave = nty + (iwave - 1)
	    do ity = 1, nty {
		if (jwave <= nwave)
		    line[jwave] = line[jwave] + asrc[itx,ity] * spec[iwave]
		jwave = jwave - 1
	    }
	}

	# Normalize to line flux

	total = asumr (line, nwave)
	if (total <= 0.0)
	    total = 1.0

	call amulkr (line, lineflux / total, line, nwave)
	return (YES)
end

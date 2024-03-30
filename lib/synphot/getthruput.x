# GETTHRUPUT -- Compute total throughput and error on requested wavelength set

procedure getthruput (mxfile, ncomp, filelist, mxparam, nparam, paramlist, 
		      nwave, wave, thruput, thruerr)

int	mxfile			# i: maximum length of file names
int	ncomp			# i: number of file names
char	filelist[mxfile,ARB]	# i: list of component througput files
int	mxparam			# i: maximum number of component parameters
int	nparam[ARB]		# i: number of component parameters
real	paramlist[mxparam,ARB]	# i: component parameters
int	nwave			# i: number of wavelengths
real	wave[ARB]		# i: wavelength grid
real	thruput[ARB]		# o: total throughput
real	thruerr[ARB]		# o: throughput error
#--
bool	noerr
int	icomp, iwave
pointer	sp, compthru, comperr

errchk	evalfilt

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (compthru, nwave, TY_REAL)
	call salloc (comperr, nwave, TY_REAL)

	# Initialize output variables

	noerr = false
	call amovkr (1.0, thruput, nwave)
	call amovkr (0.0, thruerr, nwave)

	do icomp = 1, ncomp {
	    # Compute component throughput and error on wavelength grid

	    call evalfilt (filelist[1,icomp], nparam[icomp], 
			   paramlist[1,icomp], nwave, wave, noerr, 
			   Memr[compthru], Memr[comperr])

	    # Compute product of component throughputs

	    do iwave = 1, nwave {
		if (Memr[compthru+iwave-1] <= 0) {
		    thruput[iwave] = 0.0
		} else {
		    thruput[iwave] = thruput[iwave] * Memr[compthru+iwave-1]
		}
	    }

	    # Weight component error by throughput and compute sum of squares

	    if (! noerr) {
		do iwave = 1, nwave {
		    if (Memr[compthru+iwave-1] <= 0.0) {
			thruerr[iwave] = 0.0
		    } else {
			Memr[comperr+iwave-1] = 
			    Memr[comperr+iwave-1] / Memr[compthru+iwave-1]
			thruerr[iwave] = thruerr[iwave] + 
			    Memr[comperr+iwave-1] * Memr[comperr+iwave-1] 
		    }
		}
	    }
	}

	# Multiply throughput by error for total error

	if (noerr) {
	    call amovkr (INDEFR, thruerr, nwave)
	} else {
	    do iwave = 1, nwave
		thruerr[iwave] = thruput[iwave] * sqrt (thruerr[iwave])
	}

	call sfree (sp)
end

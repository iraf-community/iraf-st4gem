# COMPCOUNTS -- Compute the counts from the spectrum and bandpass

#* HISTORY *
#* B.Simon	17-Mar-93	original

procedure compcounts (hstarea, nwave, wave, filt, spec, form, counts)

real	hstarea		# i: telescope area in cm^2
int	nwave		# i: Length of wavelength array
real	wave[ARB]	# i: wavelength array
real	filt[ARB]	# i: bandpass as a function of wavelength
real	spec[ARB]	# i: spectrum as a function of wavelength
char	form[ARB]	# i: units of spectrum
pointer	counts		# o: counts as a function of wavelength
#--
bool	status
int	i

string	badform  "compcounts: spectrum has illegal form"

int	strsearch()

begin
	# Allocate memory for counts array
	# Memory must be deallocated by calling program

	call malloc (counts, nwave, TY_REAL)

	# Compute count rate from bandpass and spectrum

	if (strsearch (form, "mag") > 0) {
	    do i = 1, nwave {
		if (filt[i] > 0.0 && ! IS_INDEFR(filt[i]) && 
		    ! IS_INDEFR(spec[i])) {
		    Memr[counts+i-1] = spec[i] - 2.5 * alog10 (filt[i])
		} else {
		    Memr[counts+i-1] = INDEFR
		}
	    }

	} else {
	    call bmulr (spec, filt, Memr[counts], nwave)
	}

	# Convert to count units

	call specform (nwave, wave, Memr[counts], form, Memr[counts], "counts",
		       status)
	if (! status)
	    call error (1, badform)

end

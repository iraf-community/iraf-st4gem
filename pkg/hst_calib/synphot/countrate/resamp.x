# RESAMP -- Resample the spectrum on a new wavelength set

#* HISTORY *
#* B.Simon	15-Apr-94	original
#* B.Simon	26-May-94	modified to call new bandwave
#* B.Simon	20-Aug-96	modified to call rebinspec

procedure resamp (wavetab, mode, cenwave, nwave, wave, band, flux)

char	wavetab[ARB]	# i: wavelength table name or coefficient set
char	mode[ARB]	# i: instrument mode string
real	cenwave		# i: central wavelength
int	nwave		# u: number of wavelengths in set
pointer	wave		# u: pointer to array of wavelengths
pointer	band		# u: pointer to passband array
pointer	flux		# u: pointer to array of fluxes
#--
int	nwave2
pointer	wave2, band2, flux2

bool	bandwave()

begin
	# Get the wavelength set appropriate for instrument and grating

	if (! bandwave (wavetab, mode, cenwave, nwave2, wave2))
	    return

	# Interpolate to resample bandpass and flux on new wavelength set

	call malloc (band2, nwave2, TY_REAL)
	call malloc (flux2, nwave2, TY_REAL)

	call syninterp (nwave, Memr[wave], Memr[band],
			nwave2, Memr[wave2], Memr[band2])

	call rebinspec (nwave, Memr[wave], Memr[flux], wavetab, 
			nwave2, Memr[wave2], Memr[flux2])

	# Free old arrays and set their pointers to the new arrays

	call mfree (band, TY_REAL)
	call mfree (wave, TY_REAL)
	call mfree (flux, TY_REAL)

	nwave = nwave2
	band = band2
	wave = wave2
	flux = flux2

end


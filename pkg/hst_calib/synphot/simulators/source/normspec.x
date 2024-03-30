#* HISTORY *
#* B.Simon	17-Feb-95	original

# NORMSPEC -- Compute normalized flux of a spectrum

procedure normspec (spectrum, exptime, mag, magband, magform, 
		    grftable, cmptable, nwave, wave, flux)

char	spectrum[ARB]	# i: Object spectrum
real	exptime		# i: exposure time
real	mag		# i: Object magnitude
char	magband[ARB]	# i: Magnitude passband
char	magform[ARB]	# i: Magnitude form
char	grftable[ARB]	# i: Instrument graph table
char	cmptable[ARB]	# i: Component name table
int	nwave		# i: Length of wavelength and thruput arrays
real	wave[ARB]	# i: Wavelength array
real	flux[ARB]	# o: Flux array
#--
int	degree
pointer	sp, band, spec, code

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (band, SZ_LINE, TY_CHAR)
	call salloc (spec, SZ_LINE, TY_CHAR)
	call salloc (code, SZ_LINE, TY_INT)

	# Create the renormalized spectrum expression

	call exp_rewrite (magband, Memc[band], SZ_LINE)

	call sprintf (Memc[spec], SZ_LINE, "%g*rn(%s,%s,%g,%s)")
	call pargr (exptime)
	call pargstr (spectrum)
	call pargstr (Memc[band])
	call pargr (mag)
	call pargstr (magform)

	# Evaluate the expression

	call syncompile (Memc[spec], Memi[code], SZ_LINE)
	call syncalc (Memi[code], SZ_LINE, NULL, nwave, wave, 
		      grftable, cmptable, flux, degree)
	
	call sfree  (sp)
end

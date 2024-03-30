# RDTHRU -- Read a component throughput from a table 

procedure rdthru (filename, nwave, wave, thruput)

char	filename[ARB]	# i: name of file containing spectrum
int     nwave           # i: length of wavelength and spectrum arrays
real    wave[ARB]       # i: wavelengths at which spectrum is computed
real	thruput[ARB]	# o: throughput at sampled wavelengths
#--
bool	noerr
int	nparam
pointer	sp, comperr
real	param[1]

errchk	evalfilt

begin
	call smark (sp)
	call salloc (comperr, nwave, TY_REAL)

	nparam = 0
	noerr = false

	call evalfilt (filename, nparam, param, nwave, wave, noerr, 
		       thruput, Memr[comperr])

	call sfree (sp)
end

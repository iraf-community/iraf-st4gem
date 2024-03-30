# EBMVFUNC -- Default interstellar reddening function

procedure ebmvfunc (extval, nwave, wave, band)

real	extval		# i: extinction value
int	nwave		# i: length of wavelength and bandpass arrays
real	wave[ARB]	# i: wavelength set output is produced on
real	band[ARB]	# o: output bandpass
#--

begin
	# Call the extended reddening function using Seaton's
	# reddening formula as a default

	call ebmvxfunc (extval, "gal1", nwave, wave, band)
end

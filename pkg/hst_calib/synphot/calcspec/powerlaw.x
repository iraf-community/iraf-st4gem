# POWERLAW -- Create a powerlaw spectrum

procedure powerlaw( script, nwv, wv, form, iw, spec )

char 	script[ARB]		# i: Command script
int	nwv			# i: number of wavelengths
real	wv[ARB]			# i: array of wavelengths
char	form[ARB]		# i: form (units) of spectrum
int	iw			# io: position in script
real 	spec[ARB]		# o: power law spectrum

char	cplwave[SZ_LINE], cplindx[SZ_LINE]
int 	ip, i, nchar
int 	ctor(), ctowrd()
real	plwave, plindx

begin

	# Reference wavelength
	nchar = ctowrd( script, iw, cplwave, SZ_LINE)
	ip = 1
	nchar = ctor( cplwave, ip, plwave)

	# Powerlaw index
	nchar = ctowrd( script, iw, cplindx, SZ_LINE)
	ip = 1
	nchar = ctor( cplindx, ip, plindx)

	if( plwave <= 0. ) {
	   call printf("Invalid powerlaw wavelength in CALCSPEC, %e.\n")
	      call pargr(plwave)
	   plwave = 1.
	}

	do i=1,nwv
	   spec[i] = ( wv[i] / plwave ) ** plindx

	call strcpy("JY",form,SZ_FNAME) 

end

# HSPEC -- Calculate hydrogen spectra

procedure hspec( script, nwv, wv, form, iw, sp)

char	script[ARB]		# i: Command script
int	nwv			# i: number of wavelengths
real 	wv[ARB]			# i: wavelength array
char 	form[ARB]		# io: form(units) of spectrum 
int	iw			# io: position in script
real 	sp[ARB]			# o: spectral array

char	ctemp[SZ_FNAME], cden[SZ_FNAME]
int	i, ip, nchar
int	ctor(), ctowrd()
real	temp, colden
real	hydnu()

begin

	# Temperature
	nchar = ctowrd( script, iw, ctemp, SZ_LINE)
	ip = 1
	nchar = ctor( ctemp, ip, temp)

	# Column density
	nchar = ctowrd( script, iw, cden, SZ_LINE)
	ip = 1
	nchar = ctor( cden, ip, colden)

	if( temp <= 0. ) {
	   call printf("** Calcspec invalid hydrogen temperature, %r\n")
	      call pargr(temp)
	}

	do i=1,nwv
	   sp[i] = hydnu( wv[i], temp, colden )

	call strcpy("JY", form, SZ_LINE)
end

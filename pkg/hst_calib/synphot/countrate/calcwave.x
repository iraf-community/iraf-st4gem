#* HISTORY *
#* B.Simon	27-Mar-96	original
#* B.Simon	18-Jul-96	calcstep extracted

# CALCWAVE -- Calculate a wavelength set from its coefficients

procedure calcwave (coef, wave, nwave)

char	coef[ARB]	# i: list of coefficients
pointer	wave		# o: wavelength set
int	nwave		# o: length of wavelength set
#--
int	iwave
real	a, b, c

begin
	# Calculate the coefficients used to generate the wavelength set

	call calcstep (coef, a, b, c, nwave)

	# Calculate the wavelength set

	call malloc (wave, nwave, TY_REAL)

	do iwave = 0, nwave-1 
	    Memr[wave+iwave] = ((a * iwave) + b) * iwave + c

end

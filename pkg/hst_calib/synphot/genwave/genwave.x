#* HISTORY *
#* E.Medeiros	01-Jan-89	original
#* B.Simon	02-Jun-94	rewritten to call waveset

# GENWAVE -- Generate a wavelength set

procedure genwave (minwave, maxwave, dwave, dvel, nwave, wave)

real	minwave		# i: minimum wavelength
real	maxwave		# i: maximum wavelength
real	dwave		# i: wavelength interval
real	dvel		# i: velocity interval
int	nwave		# o: number of wavelengths
pointer	wave		# o: wavelength set
#--
bool	logspace
real	clight, wmin, wmax, meanwave, meanvel

data	clight  / 2.997925e5 /

string	badwave  "Wavelength interval is not positive"
string	badvel   "Velocity interval is not positive"
string	wavefmt  "A/pixel = %0.3g Mean km/s/pixel = %0.3g\n"
string	velfmt   "km/s/pixel = %0.3g Mean Angstroms/pixel = %0.3g\n"
string	rngfmt   "Range = %0.4g to %0.4g Number of Pixels = %0.4d\n"

begin
	# Calculate number of points in wavelength set and recalculate
	# maximum wavelength

	if (! IS_INDEFR (dwave)) {
	    if (dwave <= 0.0)
		call printerr_real (badwave, dwave)

	    logspace = false
	    nwave = 1 + nint ((maxwave - minwave) / dwave)

	    wmin = minwave
	    wmax = minwave + dwave * (nwave - 1)

	    if (nwave <= 1)
		meanvel = 0.0
	    else
		meanvel = clight * alog (wmax / wmin) / (nwave - 1)

	    call eprintf (wavefmt)
	    call pargr (dwave)
	    call pargr (meanvel)

	} else {
	    if (dvel <= 0.0)
		call printerr_real (badvel, dvel)

	    logspace = true
	    nwave = 1 + nint (clight * alog (maxwave / minwave) / dvel)

	    wmin = minwave
	    wmax = wmin * exp ((nwave - 1) * dvel / clight)

	    if (nwave <= 1)
		meanwave = 0.0
	    else
		meanwave = (wmax - wmin) / (nwave - 1)

	    call eprintf (velfmt)
	    call pargr (dvel)
	    call pargr (meanwave)
	}

	call eprintf (rngfmt) 
	call pargr (wmin)
	call pargr (wmax)
	call pargi (nwave)

	call flush (STDOUT)

	# Calculate wavelength set from minimum, maximum wavelength
	# and number of points

	call malloc (wave, nwave, TY_REAL)
	call waveset (logspace, wmin, wmax, nwave, Memr[wave])

end

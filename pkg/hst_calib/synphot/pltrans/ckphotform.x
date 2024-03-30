# CKPHOTFORM -- Checks that data pairs have the correct forms

procedure ckphotform(nwave, wave, modes, forms, flux, err, mode, 
	             diff, rindex, nrow, valpt, npt,
	             grtbl, cmptbl, dat)

int	nwave			# i: Number of wavelengths in wave
pointer	wave			# i: Pointer to wavelength array
char	modes[SZ_LINE,ARB]	# i: Modes of flux array points
char	forms[SZ_COLUNITS,ARB]	# i: Forms of flux array points
real	flux[ARB]		# i: Flux values 
real	err[ARB]		# i: 1 sigma errors on flux
char	mode[ARB]		# i: Mode for data
bool	diff			# i: Color difference?
int	rindex[ARB]		# i: Sorted indices pointing to flux pts
int	nrow			# i: Number of points in rindex
int	valpt[2,ARB]		# i: Indices of valid points
int	npt			# i: Points in valpt
char	grtbl[ARB]		# i: Graph table
char	cmptbl[ARB]		# i: component table
real	data[2,ARB]		# o: valid data array

real	band

string	nulcmd	""

begin

	# Calculate the pivot wavelength for the first datapoint
	call synphot( nwave, wave, band, nulcmd, "pivlam",
	              modes[1,valpt[1,1]], grtbl, cmptbl, pivlam )

	dat[1,1] = flux[valpt[1,1]]
	sig[1,1] = err[valpt[1,1]]

	call sphotform(1, pivlam, dat[1,1], sig[1,1], forms[1,valpt[1,ix]],
	               dat[1,1], sig[1,1], form )

	# Check and change form of rest of x datapoints
	do ix = 2, npt {
	   dat[1,ix] = flux[valpt[1,ix]]
	   sig[1,ix] = sig[valpt[1,ix]]

	   # If form of current point is different than form of previous
	   # point then change form
	   if ( strne(forms[1,valpt[1,ix]], forms[1,validx[1,ix-1]]) )
	      call sphotform(1, pivlam, dat[1,ix], sig[1,ix], 
	                     forms[1,valid[1,ix]],
	                     dat[1,ix], sig[1,ix], form )

	}


	# Calculate the pivot wavelength for the first datapoint
	call synphot( nwave, wave, band, nulcmd, "pivlam",
	              modes[1,valid[2,1]], grtbl, cmptbl, pivlam )

	dat[2,1] = flux[valpt[2,1]]
	sig[2,1] = err[valpt[2,1]]

	call sphotform(1, pivlam, dat[2,1], sig[2,1], forms[1,valpt[2,1]],
	               dat[2,1], sig[2,1], form )

	do ix = 2, nx {
	   dat[2,ix] = flux[valpt[2,ix]]
	   sig[2,ix] = sig[valpt[2,ix]]

	   # If form of current point is different than form of previous
	   # point then change form
	   if ( strne(forms[1,valpt[2,ix]], forms[1,valpt[2,ix-1]]) )
	      call sphotform(1, pivlam, dat[2,ix], sig[2,ix], 
	                     forms[1,valpt[2,ix]],
	                     dat[1,ix], sig[2,ix], form )

	}

end

# LINTERP -- Linear interpolation, with extrapolation to constant value

#* HISTORY *
#* B.Simon	21-Sep-93	Rewritten to call synextrap

procedure linterp (form, norg, xorg, yorg, npix, xdata, ydata)

char	form[ARB]	# i: spectral form
int	norg		# i: number of original points
real	xorg[ARB]	# i: original x-values
real	yorg[ARB]	# u: original y-values
int	npix		# i: number of interpolated values
real	xdata[ARB]	# i: array of interpolated x-values
real	ydata[ARB]	# o: array of interpolated y-values
#--
int	junk
real	const

int	strsearch(), fillnull()

begin
	# Set extrapolated value to 100 for magnitude forms,
	# zero for other forms

	if (strsearch (form, "mag") > 0 || strsearch (form, "MAG") > 0) {
	    const = 100.0
	} else {
	    const = 0.0
	}

	# Remove nulls from input data

	junk = fillnull (const, norg, xorg, yorg)

	# Call interpolation routine

	call synextrap (const, norg, xorg, yorg, npix, xdata, ydata)

end

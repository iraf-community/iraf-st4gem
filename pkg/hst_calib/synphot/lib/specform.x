# SPECFORM -- Convert spectrum from one data form to another

#* HISTORY *
#* B.Simon	01-Apr-94	Rewritten to call anytophot, phottoany
#* B.Simon	26-Apr-94	Call inisyntab, clssyntab

procedure specform (nwave, wave, specin, inform, specout, outform, status)

int	nwave		# i: number of array elements
real	wave[ARB]	# i: wavelength set, in angstroms
real	specin[ARB]	# i: input flux array
char	inform[ARB]	# i: input flux form name
real	specout[ARB]	# o: output flux array
char	outform[ARB]	# i: output flux form name
bool	status		# o: was conversion successful?
#--
int	check1, check2, check3

int	is_magunit(), fillnull(), anytophot(), phottoany()
errchk	is_magunit, anytophot, phottoany

begin
	# Initialize the synphot file cache (used to read vega.dat)

	call inisyntab

	# Replace null values in input spectrum

	call amovr (specin, specout, nwave)

	if (is_magunit (inform) == NO) {
	    check1 = fillnull (0.0, nwave, wave, specout)
	} else 	{
	    check1 = fillnull (100.0, nwave, wave, specout)
	}

	# Convert form to photlam on the way to the final form

	if (check1 == OK) {
	    check2 = anytophot (inform, nwave, wave, specout)
	    check3 = phottoany (outform, nwave, wave, specout)
	}

	# Restore null values to output spectrum

	call setnull (outform, nwave, wave, specout)

	# Set status to true if both steps in conversion were performed

	status = (check2 == YES) && (check3 == YES)

	# Close the synphot file cache

	call clssyntab

end

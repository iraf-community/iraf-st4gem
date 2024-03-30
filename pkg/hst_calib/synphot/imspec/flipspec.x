#* HISTORY*
#* B.Simon	04-Mar-93	original

# FLIPSPEC -- Flip a spectrum so that the first wavelength is smallest

procedure flipspec (file, wave, spec, len)

char	file[ARB]	# i: file name (for error message)
real	wave[ARB]	# u: wavelenth array
real	spec[ARB]	# u: flux array
int	len		# i: length of arrays
#--
int	i, j, dir
pointer	sp, errmsg
real	dwave, temp

string	badorder  "Wavelengths are not monotonic (%s)"

begin
	# Allocate memory for error messages

	call smark (sp)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	# Check to make sure wavelengths are monotonic

	dir = 0
	do i = 2, len {
	    dwave = wave[i] - wave[i-1]
	    switch (dir) {
	    case -1:
		if (dwave >= 0.0) {
		    dir = 0
		    break
		}
	    case 0:
		if (dwave == 0.0) {
		    dir = 0
		    break
		} else if (dwave > 0.0) {
		    dir = 1
		} else {
		    dir = -1
		}
	    case 1:
		if (dwave <= 0.0) {
		    dir = 0
		    break
		}
	    }
	}

	# Write error message if wavelengthd are not monotonic
	# Flip the arrays if they are monotonic decreasing

	if (dir == 0) {
	    call sprintf (Memc[errmsg], SZ_LINE, badorder)
	    call pargstr (file)
	    call error (1, Memc[errmsg])

	} else if (dir == -1) {
	    i = 1
	    j = len

	    while (i < j) {
		temp = wave[i]
		wave[i] = wave[j]
		wave[j] = temp

		temp = spec[i]
		spec[i] = spec[j]
		spec[j] = temp

		i = i + 1
		j = j - 1
	    }
	}

	call sfree (sp)
end

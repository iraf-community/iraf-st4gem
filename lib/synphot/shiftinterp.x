# SHIFTINTERP -- Interpolate a new function after optionally applying a shift

procedure shiftinterp (shift, loval, hival, value, nrow, wave, 
		       lofun, hifun, fun)

bool	shift		# i: apply shift to data before interpolating?
real	loval		# i: parameter associated with low function
real	hival		# i: parameter associated with high function
real	value		# i: parameter value to interpolate at
int	nrow		# i: number of rows in interpolated functions
real	wave[ARB]	# i: dependent variable in function
real	lofun[ARB]	# i: tabulated values of low function
real	hifun[ARB]	# i: tabulated values of high function
real	fun[ARB]	# o: interpolated function
#--
pointer	lowave, hiwave, lodata, hidata

begin
	if (! shift) {
	    call funinterp (loval, hival, value, nrow, lofun, hifun, fun)

	} else if (value == loval) {
	    call amovr (lofun, fun, nrow)

	} else if (value == hival) {
	    call amovr (hifun, fun, nrow)

	} else {
	    # Allocate memory for temporary arrays

	    call malloc (lowave, nrow, TY_REAL)
	    call malloc (hiwave, nrow, TY_REAL)
	    call malloc (lodata, nrow, TY_REAL)
	    call malloc (hidata, nrow, TY_REAL)

	    if (loval == hival) {
		# Do shift but not interpolation

		call amovr (wave, Memr[hiwave], nrow)
		call aaddkr (Memr[hiwave], hival-value, Memr[hiwave], nrow)

		call syninterp (nrow, wave, hifun, nrow, Memr[hiwave], fun)

	    } else {
		# Shift function data and interpolate between shifted data

		call amovr (wave, Memr[lowave], nrow)
		call aaddkr (Memr[lowave], loval-value, Memr[lowave], nrow)

		call amovr (wave, Memr[hiwave], nrow)
		call aaddkr (Memr[hiwave], hival-value, Memr[hiwave], nrow)

		call syninterp (nrow, wave, lofun, nrow, 
				Memr[lowave], Memr[lodata])
		call syninterp (nrow, wave, hifun, nrow,
				Memr[hiwave], Memr[hidata])

		# Interpolate between arrays of shifted data

		call funinterp (loval, hival, value, nrow, Memr[lodata], 
				Memr[hidata], fun)

	    }

	    # Free temporary arrays

	    call mfree (lowave, TY_REAL)
	    call mfree (hiwave, TY_REAL)
	    call mfree (lodata, TY_REAL)
	    call mfree (hidata, TY_REAL)
	}
end

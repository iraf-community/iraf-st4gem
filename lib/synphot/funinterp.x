# FUNINTERP -- Interpolate a new function between two tabulated functions

procedure funinterp (loval, hival, value, nrow, lofun, hifun, fun)

real	loval		# i: parameter associated with low function
real	hival		# i: parameter associated with high function
real	value		# i: parameter value to interpolate at
int	nrow		# i: number of rows in interpolated functions
real	lofun[ARB]	# i: tabulated values of low function
real	hifun[ARB]	# i: tabulated values of high function
real	fun[ARB]	# o: interpolated function
#--
int	irow
real	a, b

begin
	if (loval == hival) {
	    call amovr (lofun, fun, nrow)

	} else if (value == loval) {
	    call amovr (lofun, fun, nrow)

	} else if (value == hival) {
	    call amovr (hifun, fun, nrow)

	} else {
	    a = (hival - value) / (hival - loval)
	    b = 1.0 - a

	    do irow = 1, nrow {
		if (IS_INDEFR(lofun[irow]) || IS_INDEFR(hifun[irow])) {
		    fun[irow] = INDEFR
		} else {
		    fun[irow] = a * lofun[irow] + b * hifun[irow]
		}
	    }
	}
end

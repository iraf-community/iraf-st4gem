#* HISTORY *
#* B.Simon	17-Jun-94	original

# CHOPLIM -- Set values at or outside of limits to INDEF

procedure choplim (minval, maxval, input, output, len)

real	minval		# i: minimum value
real	maxval		# i: maximum value
real	input[ARB]	# i: input array
real	output[ARB]	# o: output array
int	len		# i: length of input and output
#--
int	i

begin
	call amovr (input, output, len)

	# Set values <= minimum value to INDEF

	if (! IS_INDEFR (minval)) {
	    do i = 1, len {
		if (IS_INDEFR (output[i]))
		    next

		if (output[i] <= minval)
		    output[i] = INDEFR
	    }
	}

	# Set values >= maximum value to INDEF

	if (! IS_INDEFR (maxval)) {
	    do i = 1, len {
		if (IS_INDEFR (output[i]))
		    next

		if (output[i] >= maxval)
		    output[i] = INDEFR
	    }
	}

end

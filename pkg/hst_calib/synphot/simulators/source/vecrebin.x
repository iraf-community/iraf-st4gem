# VECREBIN -- Rebin a vector on a different grid

procedure vecrebin (niv, ivec, nov, ovec)

int	niv		# i: length of input vector 
real	ivec[ARB]	# i: input vector
int	nov		# i: length of output vector
real	ovec[ARB]	# o: output vector
#--
int	iiv, iov
pointer	sp, ipos, opos
real	step, itotal, ototal

real	asumr()

begin
	call smark (sp)
	call salloc (ipos, niv, TY_REAL)
	call salloc (opos, nov, TY_REAL)

	# Compute the positional spacing of the input and output arrays

	do iiv = 0, niv-1
	    Memr[ipos+iiv] = iiv

	step = real (niv - 1) / real (nov - 1)

	do iov = 0, nov-1
	    Memr[opos+iov] = iov * step

	# Interpolate on the new spacing

	call syninterp (niv, Memr[ipos], ivec, nov, Memr[opos], ovec)

	# Normalize the new output to the same sum as the original
	# If the new sum is zero, set the midpoint to the original sum

	itotal = asumr (ivec, niv)
	ototal = asumr (ovec, nov)

	if (ototal > 0.0) {
	    call amulkr (ovec, itotal/ototal, ovec, nov)

	} else {
	    call aclrr (ovec, nov)
	    iov = (nov + 1) / 2
	    ovec[iov] = itotal
	}

	call sfree (sp)
end


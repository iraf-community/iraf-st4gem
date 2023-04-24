# NL_ZEVAL --  Evaluates non-linear function at a specified point.
# Functions are evaluated with unscaled coefficients.

real procedure nl_zeval (nl, x, y)

pointer	nl		# i: Curve descriptor.
real	x		# i: Independent x variable.
real	y		# i: Independent y variable.

#--
real	xa[1]		# Fake arrays.
real	ya[1]
real	za[1]

errchk	nl_vector

begin
	xa[1] = x
	ya[1] = y
	call nl_vector (nl, xa, ya, za, 1)
	return (za[1])
end
                          

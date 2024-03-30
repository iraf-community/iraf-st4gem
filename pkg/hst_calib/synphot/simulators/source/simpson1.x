#* HISTORY *
#* B.Simon	11-Jul-95	Derived from simpson2

# SIMPSON1 -- Simpson integration in one dimension

procedure simpson1 (xlo, xhi, func, term, val)

real	xlo	# i: lower x boundary
real	xhi	# i: upper x boundary
real	func	# i: function to be evaluated
real	term	# u: partial term
real	val	# o: integrated value
#--
real	xmid
extern	func

begin
	# No, not O.J. Simpson !

	xmid = 0.5 * (xlo + xhi)

	# Keeping track of the partial term saves a function evaluation
	# since the last value of xhi is the current xlo

	if (IS_INDEFR (term))
	    term = func (xlo)

	val = 4. * func (xmid)
	val = val + term

	term = func (xhi)
	val = val + term

	val = val * (xhi - xlo) / 6.
end


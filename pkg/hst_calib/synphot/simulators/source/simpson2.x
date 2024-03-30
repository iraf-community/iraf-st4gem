#* HISTORY *
#* B.Simon	17-Feb-95	original

# SIMPSON2 -- Simpson integration in two dimensions

procedure simpson2 (xlo, xhi, ylo, yhi, func, sum, val)

real	xlo	# i: lower x boundary
real	xhi	# i: upper x boundary
real	ylo	# i: lower y boundary
real	yhi	# i: upper y boundary
real	func	# i: function to be evaluated
real	sum	# u: partial sum
real	val	# o: integrated value
#--
real	xmid, ymid
extern	func

begin
	# No, not O.J. Simpson !

	xmid = 0.5 * (xlo + xhi)
	ymid = 0.5 * (ylo + yhi)

	# Keeping track of the partial sum saves three function evaluations
	# since the last value of xhi is the current xlo

	if (IS_INDEFR (sum))
	    sum = func (xlo, ylo) + 4. * func (xlo, ymid) + func (xlo, yhi)

	val = 4. * func (xmid, ylo) + 16. * func (xmid, ymid) +
	      4. * func (xmid, yhi)
	val = val + sum

	sum = func (xhi, ylo) + 4. * func (xhi, ymid) + func (xhi, yhi)
	val = val + sum

	val = val * (xhi - xlo) * (yhi - ylo) / 36.
end


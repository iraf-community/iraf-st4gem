define	RANGE	5.

# NL_NGAUSS -- Returns a Gaussian-distributed number with dispersion sigma.

real procedure nl_ngauss (sigma, seed)

long	seed			# io: random number seed
real	sigma			# i : sigma of Gaussian distribution

#--
real	x1, x2, ee, ss, aa, ab

real	urand()

begin
	x2 = 1.
	ee = 0.
	ss = sqrt(2.) * sigma
	aa = RANGE * sigma
	ab = 2. * aa
	while (x2 > ee) {
	    x1 = ab * urand (seed) - aa
	    ee = exp (-((x1 / ss) ** 2))
	    x2 = urand (seed)
	}
	return (x1)
end
                                

   
# LEGENDRE -- Calculate Legendre polynomial of order n at point x

real procedure legendre( x, n )

real	x	# i: argument of polynomial
int	n	# i: order of polynomial

int	jc
real	pnm2, pnm1, pn

begin

	if ( n == 0 )
	   return (1.)

	else if (n == 1 )
	   return (x)

	# Initialize P0 and P1
	pnm2 = 1.
	pnm1 = x

	for (jc = 2; jc <= n; jc=jc+1) {

	   pn = ( (2*jc - 1) * x * pnm1 - (n-1) * pnm2 )/jc
	   pnm2 = pnm1
	   pnm1 = pn

	}

	return ( pn )

end

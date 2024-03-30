# BNU -- Planck function
# input:
#	wave	r4 wavelength in angstroms
#	temp	r4 kelvin temperature
# output:
#	bnu	r4 blackbody intensity (erg/cm2/s/hz/ster)
#--
# jun 1989 Heith Horne @ stsci - clean up old subroutine
# jul 1989 Dave Bazell -- SPP version

real procedure bnu( wave, temp )

real wave		# input wavelength in Angstroms
real temp		# input temperature in Kelvin
real x			# internal variable that changes, related to wave*temp
real factor		# Bose factor; number of modes with wave at temp

real c1			# hc/k [Angstroms-Kelvin]
real c2			# numerical factor to make units come out OK
data c1/1.43883e8/
data c2/1.95722e5/

begin

	x = wave * temp
	if( x <= 0. ) return(0.)
	x = c1 / x
        if( x < 1.e-4 )
	   factor = 2. / ( x * ( x + 2. ) )
        else if( x < 85. )
	   factor = 1. / ( exp( x ) - 1. )
        else
	   return(0.)

	x = x * temp / c2
	return( x * x * x * factor )

end 

# HYDNU -- LTE hydrogen slab
#
#  emission from hydrogen.  bound-free and free-free 
#
# Input:
#	WAVE	= Wavelength (Angstroms)
#	TEMP	= Kelvin Temperature
#	COLDEN	= column density (cm^-3) or log(column density)
# Output:
#	HYDNU	= LTE source function (erg/s/cm^2/Hz/ster)
#--
# Nov 1985 KDH @ STScI
# Mar 1989 KDH @ STScI - use column density
# Jul 1989 SPP version by Dave Bazell

real procedure hydnu( wave, temp, colden )

real	bnu()		#
real	ophyd()		#
real	stimcor()	# 
real	wave		# wavelenth for calculation
real	temp		# temperature of gas
real	opacity		# calculated opacity of slab
real	colden		# calculated column density
real	tau		# calculated optical depth
real	factor
real	z		#
real	value		# function value to be returned

begin
	value=0.
	if ( wave <= 0. || temp <= 0. ) return

# source function
	value = bnu( wave, temp)
	if (value <= 0. ) return

# opacity
	z = 1.
	opacity = ophyd( wave, temp, z ) * stimcor( wave, temp )
	if ( opacity <= 0. ) return

# optical depth
	if ( colden > 0. && colden < 80. )
	   tau = opacity * 10**colden
	else
	   tau = opacity * colden

# radiative transfer
	if ( tau < 1.e-3 )
	   factor = tau * ( 1. - tau*(0.5 - tau/6.) )
	else
	   factor = 1. - exp(-tau)
	return (value * factor)
end

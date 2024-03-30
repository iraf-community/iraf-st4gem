include "constants.h"

# STIMCOR -- Evaluates stimulated emission factor 1 - exp(-hnu/kt)
# input:
#	wave	r4 wavelength (angstroms)
#	temp	r4 kelvin temperature
# output:
#	stimcor	r4 (1 - exp(-hnu/kt))
#--
# may 1989 kdh @ stsci - revise comments
# jul 1989 kdh @ stsci - use cgs to get physical constants
# jul 1989 Dave Bazell -- SPP version

define FACTOR HPLANCK * CLIGHT * 1.e8 / KBOLTZ

real procedure stimcor( wave, temp )

real	wave			# wavelength in Angstroms 
real	temp			# temperature in Kelvin

begin
	if( wave <= 0. || temp <= 0. )
	   return(1.)
	else
	   return (1. - exp( -FACTOR/(temp*wave) ) )
end

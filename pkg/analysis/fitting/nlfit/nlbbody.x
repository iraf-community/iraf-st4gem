include <mach.h>
include "nlfit.h"

define	ARGMIN	      1.0D-3	# Threshold for long wavelength approximation
define	ARGEMAX	      88.D+0	# Maximum argument of exp()

# NL_BBODY -- planck function evaluation.
# Given the temperature, this routine computes the blackbody 
# intensity at each data point. Double precision is used throughout.  
# The formula used is:
#      
#                        ref**5 * (Exp(efactor/ref) - 1)
#       ff(i) = bbref * ________________________________
#                       x[i]**5 * (Exp(efactor/x[i] - 1)
#
#   where  efactor = c1 * xscale / temp
#
#   and x[i] and ref are input directly in the formula in case of wavelength 
#   units, or as x[i] = 1./ x[i] and ref = 1./ ref in case of frequency or
#   energy units.
#
#   The constant c1 depends on the physical units of the x variable. The
#   following are supported:
#
#        x in cm:	  c1 = 1.438786D+0
#        x in meter:	  c1 = 1.438786D-2
#        x in Angstrom:   c1 = 1.438786D+8
#        x in Hz:         c1 = 4.795953D-11
#        x in keV         c1 = 1.159650D+7
#      
#   X scaling:
#   ________________
#   The x array in the calling program may be scaled to arbitrary units of 
#   wavelength or frequency. Because physical units appear explicitly 
#   in the Planck function formula, it is necessary to convert the scaled
#   x values (and ref) back to unscaled ones. Argument xscale, which is the 
#   number of scaled units in one physical unit, is required for this purpose:
#
#          physical x = scaled x / xscale
#

procedure nl_bbody (x, z, ndata, temp, bbref, ref, xscale, unit)

real	x[ARB]		# i: Independent variable * xscale
real	z[ARB]		# o: Black-body spectrum (joules/meter**2/sec/meter)
int	ndata		# i: Number of data points.
real	temp		# i: Temperature in Kelvin
real	bbref		# i: Brightness at reference x
double	ref		# i: Reference x
double	xscale		# i: X scaling factor
int	unit		# i: X physical unit

#--

double	arg, argbuf
double	efactor, numer, denon
double	c1, ref1, xsc1, x1
int	i, iunit

begin
	iunit = unit
	if (temp <= EPSILON) {
	    call aclrr (z, ndata)
	} else {
	    ref1 = ref
	    xsc1 = xscale
	    switch (unit) {
	        case CM:
	            c1 = 1.438786D+0
	        case METER:
	            c1 = 1.438786D-2
	        case ANGSTROM:
	            c1 = 1.438786D+8
	        case HZ:
	            c1 = 4.795953D-11
	            ref1 = 1.0D0 / ref
	            xsc1 = 1.0D0 / xscale
	        case KEV:
	            c1 = 1.15965D+7
	            ref1 = 1.0D0 / ref
	            xsc1 = 1.0D0 / xscale
	            iunit = HZ
	        default:
	            call error (0, "Non-specified physical units")
	    }

	    efactor = c1 * xsc1 / temp
	    if ((efactor/ref1) > ARGEMAX)
	        call error (0, "Exp overflow in Planck function.")
	    numer = bbref * ref1 ** 5 * (exp(efactor/ref1) - 1.0D0)

	    do i = 1, ndata {
	        if (x[i] > 0.) {
	            x1 = x[i]
	            if (iunit == HZ) 
	                x1 = 1.d0 / x1
	            arg = efactor / x1
	            denon = x1 ** 5
	            if (arg < ARGMIN) 
	                denon = denon * arg * (1.D0 + 0.5D0*arg)
	            else {
	                argbuf = min(arg, ARGEMAX)
	                denon = denon * (exp(argbuf) - 1.D0)
	            }
	            z[i] = numer / denon
	            if (arg > ARGEMAX)
	                z[i] = 0.
	        } else
	            z[i] = 0.
	    }
	}
end
                                         
                                                               
                                                             

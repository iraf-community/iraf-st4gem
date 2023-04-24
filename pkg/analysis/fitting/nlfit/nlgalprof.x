include	"nlfit.h"

define		LARGE_ARG	80.d0		# largest exp argument

# NL_GALBULGE -- Evaluates galaxy bulge brightness profile.
# A deVaucouleurs law is used. A background value is added to the profile.
# Input is asumed to be distance in either linear or r**1/4 units, 
# output is linear intensity.

procedure nl_galbulge (x, z, ndata, se, re, backgr, unit)

real	x[ARB]		# i: Indep. variable array.
real	z[ARB]		# o: Computed user function.
int	ndata		# i: Number of data points.
real	se		# i: bulge brigthness at effective radius
real	re		# i: bulge effective radius
real	backgr		# i: background level
int	unit		# i: units of distance variable

#--
int	i
double	xx, earg

begin
	do i = 1, ndata {
	    if (unit == FOURTH_ROOT)
	        xx = double (x[i]) ** 4.
	    else
	        xx = double (x[i])
	    if (re > 0.)
	        earg = -7.668d0 * ((xx / re)**0.25 - 1.)
	    else
	        earg = 0.d0
	    if (earg > LARGE_ARG)
	        earg = LARGE_ARG
	    if (earg < -LARGE_ARG)
	        earg = -LARGE_ARG
	    z[i] = se * exp(earg) + backgr
	}
end


# NL_GALDISK -- Evaluates galaxy disk brightness profile.
# The profile is an exponential disk whith hole. A background value is 
# added to the profile. Input is asumed to be distance in either linear 
# or r**1/4 units, output is linear intensity.

# The sign of (r1/xx)**3 was changed to minus by PEH, 12 Aug 1992.

procedure nl_galdisk (x, z, ndata, s0, r0, r1, backgr, unit)

real	x[ARB]		# i: Indep. variable array.
real	z[ARB]		# o: Computed user function.
int	ndata		# i: Number of data points.
real	s0		# i: disk central brigthness
real	r0		# i: disk scale
real	r1		# i: disk "hole" radius
real	backgr		# i: backgr level
int	unit		# i: units of distance variable

#--
int	i
double	xx, earg

begin
	do i = 1, ndata {
	    if (unit == FOURTH_ROOT)
	        xx = double (x[i]) ** 4.
	    else
	        xx = double (x[i])
	    if ((xx > 0.d0) && (r0 > 0.))
	        earg = -(double (r1) / xx)**3 - xx / double (r0)
	    else
	        earg = -LARGE_ARG
	    if (earg > LARGE_ARG)
	        earg = LARGE_ARG
	    if (earg < -LARGE_ARG)
	        earg = -LARGE_ARG
	    z[i] = s0 * exp(earg) + backgr
	}
end



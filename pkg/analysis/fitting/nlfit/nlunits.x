include "nlfit.h"

# NL_UNITS -- Set physical units.
#
# In case of Planck, twobbody or composite function, examine input
# vector and x reference, and set physical units control variable.
# In case of galaxy brigthness profile, examine input vector and
# set in linear or r**1/4 mode. Also check for incompatibilities.

int procedure nl_units (nl, x, n)

pointer nl			# i: Curve descriptor.
real	x[ARB]			# i: Input vector.
int	n			# i: Last datum on vector.

#--
int	i, j
int	nl_xunit()
real	maxx

errchk	nl_xunit

begin
	if ((NL_FITFUNC(nl) == BBODY) || (NL_FITFUNC(nl) == COMPOSITE) ||
	    (NL_FITFUNC(nl) == TWOBBODY)) {
	    if (NL_UNIT(nl) == AUTO) {
	        i = nl_xunit (abs(x[n/2]))
	        if (NL_FITFUNC(nl) == COMPOSITE)
	            j = nl_xunit (abs(Memr[NL_PARAMS(nl)+5]))
	        else
	            j = nl_xunit (abs(Memr[NL_PARAMS(nl)+2]))
	        if ((i == ERR) || (j == ERR)) {
	            call eprintf ("Invalid physical units in Planck function.\n")
	            return (ERR)
	        }
		if (i != j) {
	            call eprintf ("Data and ref. not in the same units.\n")
	            return (ERR)
	        }
	        NL_UNIT(nl) = i
	    } else if ((NL_UNIT(nl) != ANGSTROM) && (NL_UNIT(nl) != METER) &&
		       (NL_UNIT(nl) != CM)       && (NL_UNIT(nl) != HZ)    &&
		       (NL_UNIT(nl) != KEV)) {
	        call error (0, "xunit incompatible with function type.\n")
	    }
	} else if ((NL_FITFUNC(nl) == GALPROF) || (NL_FITFUNC(nl) == GALBULGE)
		|| (NL_FITFUNC(nl) == GALDISK)) {
	    if (NL_UNIT(nl) == AUTO) {
	        maxx = -1.e5
	        do i = 1, n {
	            if (x[i] > maxx)
	                maxx = x[i]
	        }
	        if (maxx < 6.7)			# this is a 2000 pixel radius
	            NL_UNIT(nl) = FOURTH_ROOT
	        else
	            NL_UNIT(nl) = LINEAR
	    } else if ((NL_UNIT(nl) != FOURTH_ROOT) && 
		       (NL_UNIT(nl) != LINEAR)) {
	        call error (0, "xunit incompatible with function type.\n")
	    }
	}
	return (OK)
end
                                           
           

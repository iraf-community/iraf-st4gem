include "nlfit.h"

#  NL_DEFC -- Returns appropriate coefficient value for computing function
#             derivative when the input coefficient is exactly zero.
#
#  The derivative computation algorithm in nl_fdev breaks down when the
#  coefficient value is exactly zero. This routine returns a value close
#  to zero, but how close depends on the particular function type.
#  This routine will not fix erroneous cases, e.g. setting the reference
#  wavelength of a black body function to zero. It merely minimizes the 
#  chance of numerical problems in the derivative computation.


real procedure nl_defc (nl, ic, x, z)

pointer nl		# i: Curve descriptor.
int	ic		# i: Coefficient index (1-indexed)
real	x		# i: Typical value of independent variable(s)
real	z		# i: Typical value of dependent variable
#--
int	of, i
real	rval

begin
	# Compute coefficient offset in coefficient array.
	of = ic - 1

	switch (NL_FITFUNC(nl)) {

	case BBODY: 
	    if (of == NL_BTEMP) rval = 1000.
	    if (of == NL_BAMPL) rval = z / 1000.
	    if (of == NL_BREF)  rval = x / 1000.

	case POWERLAW:
	    if (of == NL_PINDEX) rval = 1.E-3
	    if (of == NL_PAMPL)  rval = z / 1000.
	    if (of == NL_PREF)   rval = x / 1000.

	case COMPOSITE:
	    if (of == NL_CINDEX) rval = 1.E-3
	    if (of == NL_CTEMP)  rval = 1000.
	    if (of == NL_CPAMPL) rval = z / 1000.
	    if (of == NL_CBAMPL) rval = z / 1000.
	    if (of == NL_CPREF)  rval = x / 1000.
	    if (of == NL_CBREF)  rval = x / 1000.

	case TWOBBODY:
	    if (of == NL_TTEMP1) rval = 1000.
	    if (of == NL_TTEMP2) rval = 1000.
	    if (of == NL_TAMPL1) rval = z / 1000.
	    if (of == NL_TAMPL2) rval = z / 1000.
	    if (of == NL_TREF1)  rval = x / 1000.
	    if (of == NL_TREF2)  rval = x / 1000.

	case GALBULGE:
	    if (of == NL_GBME)    rval = z / 1000.
	    if (of == NL_GBRE)    rval = 0.5
	    if (of == NL_GBACKGR) rval = z / 1000.

	case GALDISK:
	    if (of == NL_GDM0)    rval = z / 1000.
	    if (of == NL_GDR0)    rval = 0.5
	    if (of == NL_GDRH)    rval = 0.5
	    if (of == NL_GBACKGR) rval = z / 1000.

	case GALPROF:
	    if (of == NL_GDM0)    rval = z / 1000.
	    if (of == NL_GDR0)    rval = 0.5
	    if (of == NL_GDRH)    rval = 0.5
	    if (of == NL_GBACKGR) rval = z / 1000.
	    if (of == NL_GBME)    rval = z / 1000.
	    if (of == NL_GBRE)    rval = 0.5
	    if (of == NL_GBACKGR) rval = z / 1000.

	case GAUSSIANS:
	    if (of == NL_GA) rval = z / 1000.
	    if (of == NL_GB) rval = z / 1.E6
	    do i = 1, NL_NPAR(nl) / 3 {
	        if (of == NL_GAMPL(i))  rval = z / 1000.
	        if (of == NL_GCENT(i))  rval = x / 1000.
	        if (of == NL_GFWHM(i)) rval = x / 1000.
	    }

	case CGAUSS:
	    if (of == NL_GA)       rval = z / 1000.
	    if (of == NL_GB)       rval = z / 1.E6
	    if (of == NL_GAMPL(1)) rval = z / 1000.
	    if (of == NL_GCENT(1)) rval = x / 1000.
	    if (of == NL_GFWHM(1)) rval = x / 1000.
	    if (NL_NPAR(nl) > 5) {
	        do i = 2, NL_NPAR(nl) / 3 {
	            if (of == NL_GAMPL(i)) rval = 1.E-5
	            if (of == NL_GCENT(i)) rval = x / 1.E6
	            if (of == NL_GFWHM(i)) rval = x / 1.E6
	        }
	    }

	case USER:
	    rval = 0.

	case TWODGAUSS:
	    if (of == NL_G2A)    rval = z / 1000.
	    if (of == NL_G2AMPL) rval = z / 1000.
	    if (of == NL_G2XC)   rval = x / 1000.
	    if (of == NL_G2YC)   rval = x / 1000.
	    if (of == NL_G2ELL)  rval = 1.E-4
	    if (of == NL_G2TETA) rval = 1.E-4
	    if (of == NL_G2FWHM) rval = x / 1000.

	default:
	    call error (0, "Unknown function.")
	}

	return (rval)
end

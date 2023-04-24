include "nlfit.h"

# NL_STARTV --  Set starting values for non-defined coefficients.
#
# If there are non-defined (INDEF) coefficients, this procedure tries
# to figure out sensible starting values. They are computed by examining
# the range of independent and dependent variables.
#
# In the case of GAUSSIANS or CGAUSS functions, this routine assumes that
# the coefficient array was allocated for storage of 1 Gaussian.
#
# No guesses can be made for the USER function type, because there are no 
# clues about the meaning of each parameter.

procedure nl_startv (nl, xfit, yfit, zfit)

pointer nl			# i:  Curve descriptor.
real	xfit[ARB]		# i:  Independant x variable array.
real	yfit[ARB]		# i:  Independant y variable array.
real	zfit[ARB]		# i:  Dependant variable array.

#--
pointer	coeff
int	i, npar, npts, dim
real	xcenter, xrange, ycenter, zpeak, zbase, xpk, ypk

int	nl_stati()

errchk	nl_bbody, nl_power, nl_1gauss

begin
	# Get data from curve fitting descriptor.
	npts = nl_stati (nl, "npts")
	npar = nl_stati (nl, "npar")
	call malloc (coeff, npar, TY_REAL)
	call nl_gcoeff (nl, Memr[coeff], npar)
	dim = nl_stati (nl, "dimension")

	# Compute approximate center of X and Y distribution. 
	# This assumes data is sorted and both extremes are 
	# not INDEF.
	xrange = xfit[npts] - xfit[1]
	xcenter = 0.5 * xrange + xfit[1]
	if (dim == 2)
	    ycenter = 0.5 * (yfit[npts] - yfit[1]) + yfit[1]
	else
	    ycenter = 0.0

	# Locate max and min of independent variable.
	zpeak = -1.0E+30
	zbase = -zpeak
	xpk   = xfit[npts/2]
	ypk   = yfit[npts/2]
	do i = 1, npts {
	    if (!IS_INDEFR (zfit[i])) {
	        if (zfit[i] > zpeak) {
	            zpeak = zfit[i]
	            xpk = xfit[i]
	            if (dim == 2)
	                ypk = yfit[i]
	        }
	        if (zfit[i] < zbase)
	            zbase = zfit[i]
	    }
	}

	# Handle each case.
	switch (nl_stati (nl, "fitfunc")) {

	case POWERLAW: 
	    if (IS_INDEFR (Memr[coeff+NL_PINDEX]))
	        Memr[coeff+NL_PINDEX] = 1.0
	    if (IS_INDEFR (Memr[coeff+NL_PAMPL]))
	        Memr[coeff+NL_PAMPL] = zpeak
	    if (IS_INDEFR (Memr[coeff+NL_PREF]))
	        Memr[coeff+NL_PREF] = xcenter

	case BBODY: 
	    if (IS_INDEFR (Memr[coeff+NL_BTEMP]))
	        Memr[coeff+NL_BTEMP] = 5000.0
	    if (IS_INDEFR (Memr[coeff+NL_BAMPL]))
	        Memr[coeff+NL_BAMPL] = zpeak
	    if (IS_INDEFR (Memr[coeff+NL_BREF]))
	        Memr[coeff+NL_BREF] = xcenter

	case COMPOSITE:
	    if (IS_INDEFR (Memr[coeff+NL_CTEMP]))
	        Memr[coeff+NL_CTEMP] = 5000.0
	    if (IS_INDEFR (Memr[coeff+NL_CBAMPL]))
	        Memr[coeff+NL_CBAMPL] = zpeak / 2.0
	    if (IS_INDEFR (Memr[coeff+NL_CBREF]))
	        Memr[coeff+NL_CBREF] = xcenter
	    if (IS_INDEFR (Memr[coeff+NL_CINDEX]))
	        Memr[coeff+NL_CINDEX] = 1.0
	    if (IS_INDEFR (Memr[coeff+NL_CPAMPL]))
	        Memr[coeff+NL_CPAMPL] = zpeak
	    if (IS_INDEFR (Memr[coeff+NL_CPREF]))
	        Memr[coeff+NL_CPREF] = xcenter

	case TWOBBODY:
	    if (IS_INDEFR (Memr[coeff+NL_TTEMP2]))
	        Memr[coeff+NL_TTEMP2] = 3000.0
	    if (IS_INDEFR (Memr[coeff+NL_TAMPL2]))
	        Memr[coeff+NL_TAMPL2] = zpeak / 2.0
	    if (IS_INDEFR (Memr[coeff+NL_TREF2]))
	        Memr[coeff+NL_TREF2] = xcenter
	    if (IS_INDEFR (Memr[coeff+NL_TTEMP1]))
	        Memr[coeff+NL_TTEMP1] = 6000.0
	    if (IS_INDEFR (Memr[coeff+NL_TAMPL1]))
	        Memr[coeff+NL_TAMPL1] = zpeak
	    if (IS_INDEFR (Memr[coeff+NL_TREF1]))
	        Memr[coeff+NL_TREF1] = xcenter

	case GALBULGE:
	    if (IS_INDEFR (Memr[coeff+NL_GBRE]))
	        Memr[coeff+NL_GBRE] = xrange / 50.0
	    if (IS_INDEFR (Memr[coeff+NL_GBME]))
	        Memr[coeff+NL_GBME] = zpeak
	    if (IS_INDEFR (Memr[coeff+NL_GBACKGR]))
	        Memr[coeff+NL_GBACKGR] = zbase

	case GALDISK:
	    if (IS_INDEFR (Memr[coeff+NL_GDM0]))
	        Memr[coeff+NL_GDM0] = zpeak
	    if (IS_INDEFR (Memr[coeff+NL_GDR0]))
	        Memr[coeff+NL_GDR0] = xrange / 50.0
	    if (IS_INDEFR (Memr[coeff+NL_GDRH]))
	        Memr[coeff+NL_GDRH] = 0.0
	    if (IS_INDEFR (Memr[coeff+NL_GBACKGR]))
	        Memr[coeff+NL_GBACKGR] = zbase

	case GALPROF:
	    if (IS_INDEFR (Memr[coeff+NL_GBRE]))
	        Memr[coeff+NL_GBRE] = xrange / 50.0
	    if (IS_INDEFR (Memr[coeff+NL_GBME]))
	        Memr[coeff+NL_GBME] = zpeak
	    if (IS_INDEFR (Memr[coeff+NL_GDM0]))
	        Memr[coeff+NL_GDM0] = zpeak
	    if (IS_INDEFR (Memr[coeff+NL_GDR0]))
	        Memr[coeff+NL_GDR0] = xrange / 50.0
	    if (IS_INDEFR (Memr[coeff+NL_GDRH]))
	        Memr[coeff+NL_GDRH] = 0.0
	    if (IS_INDEFR (Memr[coeff+NL_GBACKGR]))
	        Memr[coeff+NL_GBACKGR] = zbase

	case GAUSSIANS,CGAUSS:
	    if (IS_INDEFR (Memr[coeff+NL_GB]))
	        Memr[coeff+NL_GB] = 0.0
	    if (IS_INDEFR (Memr[coeff+NL_GA]))
	        Memr[coeff+NL_GA] = zbase
	    if (IS_INDEFR (Memr[coeff+NL_GAMPL(1)]))
	        Memr[coeff+NL_GAMPL(1)] = zpeak
	    if (IS_INDEFR (Memr[coeff+NL_GCENT(1)]))
	        Memr[coeff+NL_GCENT(1)] = xpk
	    if (IS_INDEFR (Memr[coeff+NL_GFWHM(1)]))
	        Memr[coeff+NL_GFWHM(1)] = xrange / 100.

	case TWODGAUSS:
	    if (IS_INDEFR (Memr[coeff+NL_G2A]))
	        Memr[coeff+NL_G2A] = zbase
	    if (IS_INDEFR (Memr[coeff+NL_G2AMPL]))
	        Memr[coeff+NL_G2AMPL] = zpeak
	    if (IS_INDEFR (Memr[coeff+NL_G2XC]))
	        Memr[coeff+NL_G2XC] = xpk
	    if (IS_INDEFR (Memr[coeff+NL_G2YC]))
	        Memr[coeff+NL_G2YC] = ypk
	    if (IS_INDEFR (Memr[coeff+NL_G2FWHM]))
	        Memr[coeff+NL_G2FWHM] = xrange / 100.
	    if (IS_INDEFR (Memr[coeff+NL_G2ELL]))
	        Memr[coeff+NL_G2ELL] = 0.1
	    if (IS_INDEFR (Memr[coeff+NL_G2TETA]))
	        Memr[coeff+NL_G2TETA] = 20.0

	case USER:
	    do i = 1, npar {
	        if (IS_INDEFR (Memr[coeff+i-1]))
	            call error (0, 
                    "Cannot guess coefficients for this function.")
	    }

	}

	# Store guess coefficients back into curve fitting descriptor.
	call nl_scoeff (nl, Memr[coeff], npar)
	call mfree (coeff, TY_REAL)
end
        
                                                               
                     

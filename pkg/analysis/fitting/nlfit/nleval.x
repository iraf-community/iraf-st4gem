include "nlfit.h"

# NL_EVAL --  Non-linear function evaluation.
#
# This procedure evaluates a non-linear function for an array of values of 
# the independent variable. Scaled coefficients must be passed in array coeff.

procedure nl_eval (nl, x1, y1, ff, coeff, npts)

pointer	nl		# i: Curve descriptor.
real	x1[ARB]		# i: Independent x variable array.
real	y1[ARB]		# i: Independent y variable array.
real	ff[ARB]		# o: Array containing the evaluated function.
real	coeff[ARB]	# i: Coefficients.
int	npts		# i: Number of data points

#--
pointer	sp, aux
int	i
real	center, ampl

errchk	nl_bbody, nl_power, nl_1gauss, nl_userf

begin
	switch (NL_FITFUNC(nl)) {

	case BBODY: 
	    call nl_bbody (x1, ff, npts, coeff[NL_BTEMP+1], 
					 coeff[NL_BAMPL+1], 
					 double(coeff[NL_BREF+1]), 
			   double(NL_XSCALE(nl)), NL_UNIT(nl))

	case POWERLAW: 
	    call nl_power (x1, ff, npts, coeff[NL_PINDEX+1], 
					 coeff[NL_PAMPL+1], 
					 double(coeff[NL_PREF+1]))

	case COMPOSITE:
	    call smark (sp)
	    call salloc (aux, npts, TY_REAL)
	    call nl_bbody (x1, Memr[aux], npts, coeff[NL_CTEMP+1], 
						coeff[NL_CBAMPL+1], 
						double(coeff[NL_CBREF+1]), 
			   double(NL_XSCALE(nl)), NL_UNIT(nl))
	    call nl_power (x1, ff, npts, coeff[NL_CINDEX+1], 
					 coeff[NL_CPAMPL+1], 
					 double(coeff[NL_CPREF+1]))
	    call aaddr (ff, Memr[aux], ff, npts)
	    call sfree (sp)

	case TWOBBODY:
	    call smark (sp)
	    call salloc (aux, npts, TY_REAL)
	    call nl_bbody (x1, Memr[aux], npts, coeff[NL_TTEMP2+1], 
						coeff[NL_TAMPL2+1], 
						double(coeff[NL_TREF2+1]), 
			   double(NL_XSCALE(nl)), NL_UNIT(nl))
	   call nl_bbody (x1, ff, npts, coeff[NL_TTEMP1+1], 
					coeff[NL_TAMPL1+1], 
					double(coeff[NL_TREF1+1]),
			   double(NL_XSCALE(nl)), NL_UNIT(nl))
	    call aaddr (ff, Memr[aux], ff, npts)
	    call sfree (sp)

	case GALBULGE:
	    call nl_galbulge (x1, ff, npts, coeff[NL_GBME+1], coeff[NL_GBRE+1], 
			      coeff[NL_GBACKGR+1], NL_UNIT(nl))

	case GALDISK:
	    call nl_galdisk (x1, ff, npts, coeff[NL_GDM0+1], coeff[NL_GDR0+1],
			     coeff[NL_GDRH+1], coeff[NL_GBACKGR+1], NL_UNIT(nl))

	case GALPROF:
	    call smark (sp)
	    call salloc (aux, npts, TY_REAL)
	    # Disk component. Half background is added here.
	    call nl_galdisk (x1, Memr[aux], npts, coeff[NL_GDM0+1], 
			     coeff[NL_GDR0+1], coeff[NL_GDRH+1], 
			     coeff[NL_GBACKGR+1]/2., NL_UNIT(nl))
	    # Bulge component. Half background is added here.
	    call nl_galbulge (x1, ff, npts, coeff[NL_GBME+1], coeff[NL_GBRE+1], 
			      coeff[NL_GBACKGR+1]/2., NL_UNIT(nl))
	    call aaddr (ff, Memr[aux], ff, npts)
	    call sfree (sp)

	case GAUSSIANS:
	    call amulkr  (x1, coeff[NL_GB+1], ff, npts)		# linear
	    call aaddkr (ff, coeff[NL_GA+1], ff, npts)		# baseline
	    do i = 1, NL_NPAR(nl) / 3
	        call nl_1gauss (x1, ff, npts, coeff[NL_GAMPL(i)+1], 
					      coeff[NL_GCENT(i)+1],
					      coeff[NL_GFWHM(i)+1])

	case CGAUSS:
	    call amulkr  (x1, coeff[NL_GB+1], ff, npts)
	    call aaddkr (ff, coeff[NL_GA+1], ff, npts)
	    call nl_1gauss (x1, ff, npts, coeff[NL_GAMPL(1)+1], 
					  coeff[NL_GCENT(1)+1],
					  coeff[NL_GFWHM(1)+1])
	    if (NL_NPAR(nl) > 5) {
	        do i = 2, NL_NPAR(nl) / 3 {
	            ampl   = coeff[NL_GAMPL(1)+1] * coeff[NL_GAMPL(i)+1]
	            center = coeff[NL_GCENT(1)+1] + coeff[NL_GCENT(i)+1]
	            call nl_1gauss (x1, ff, npts, ampl, center, 
				    coeff[NL_GFWHM(i)+1])
	        }
	    }
	case USER:
	    call nl_userf (nl, x1, ff, Memr[NL_SY(nl)], npts, 
			   coeff, NL_NPAR(nl))

	case TWODGAUSS:
	    call nl_2dgauss (x1, y1, ff, npts, coeff[NL_G2AMPL+1], 
					       coeff[NL_G2XC+1],
					       coeff[NL_G2YC+1],
					       coeff[NL_G2FWHM+1],
					       coeff[NL_G2ELL+1],
					       coeff[NL_G2TETA+1])
	    call aaddkr (ff, coeff[NL_G2A+1], ff, npts)		# baseline

	default:
		call error (0, "Error in function specification")
	}
end
        
                                                               
                     

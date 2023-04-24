include "nlfit.h"

# NL_SCALE -- Scale function coefficients in accordance with the scaling
# performed on data values.

procedure nl_scale (nl)

pointer nl			# i: Curve descriptor.

#--
int	i

begin
	# The scaling depends on the particular type of fitting function.

	switch (NL_FITFUNC(nl)) {

	case BBODY: 
	    Memd[NL_PARSCALE(nl)+NL_BTEMP] = 1.0d0	  # Temp not scaled
	    Memd[NL_PARSCALE(nl)+NL_BAMPL] = double(NL_ZSCALE(nl))  # Amplitude
	    Memd[NL_PARSCALE(nl)+NL_BREF]  = double(NL_XSCALE(nl))  # Reference

	case POWERLAW:
	    Memd[NL_PARSCALE(nl)+NL_PINDEX] = 1.0d0	  # Index not scaled
	    Memd[NL_PARSCALE(nl)+NL_PAMPL]  = double(NL_ZSCALE(nl))  # Amplitude
	    Memd[NL_PARSCALE(nl)+NL_PREF]   = double(NL_XSCALE(nl))  # Reference

	case COMPOSITE:
	    Memd[NL_PARSCALE(nl)+NL_CINDEX] = 1.0d0	 # Index and temperature
	    Memd[NL_PARSCALE(nl)+NL_CTEMP]  = 1.0d0	 # not scaled
	    Memd[NL_PARSCALE(nl)+NL_CPAMPL] = double(NL_ZSCALE(nl)) # Amplitudes
	    Memd[NL_PARSCALE(nl)+NL_CBAMPL] = double(NL_ZSCALE(nl)) 
	    Memd[NL_PARSCALE(nl)+NL_CPREF]  = double(NL_XSCALE(nl)) # References
	    Memd[NL_PARSCALE(nl)+NL_CBREF]  = double(NL_XSCALE(nl)) 

	case TWOBBODY:
	    Memd[NL_PARSCALE(nl)+NL_TTEMP1] = 1.0d0	  # Temperatures
	    Memd[NL_PARSCALE(nl)+NL_TTEMP2] = 1.0d0	  # not scaled
	    Memd[NL_PARSCALE(nl)+NL_TAMPL1] = double(NL_ZSCALE(nl)) # Amplitudes
	    Memd[NL_PARSCALE(nl)+NL_TAMPL2] = double(NL_ZSCALE(nl)) 
	    Memd[NL_PARSCALE(nl)+NL_TREF1]  = double(NL_XSCALE(nl)) # References
	    Memd[NL_PARSCALE(nl)+NL_TREF2]  = double(NL_XSCALE(nl)) 

	case GALBULGE:
	    Memd[NL_PARSCALE(nl)+NL_GBME]    = double(NL_ZSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GBRE]    = double(NL_XSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GBACKGR] = double(NL_ZSCALE(nl))

	case GALDISK:
	    Memd[NL_PARSCALE(nl)+NL_GDM0]    = double(NL_ZSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GDR0]    = double(NL_XSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GDRH]    = double(NL_XSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GBACKGR] = double(NL_ZSCALE(nl))

	case GALPROF:
	    Memd[NL_PARSCALE(nl)+NL_GDM0]    = double(NL_ZSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GDR0]    = double(NL_XSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GDRH]    = double(NL_XSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GBME]    = double(NL_ZSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GBRE]    = double(NL_XSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GBACKGR] = double(NL_ZSCALE(nl))

	case GAUSSIANS:
	    Memd[NL_PARSCALE(nl)+NL_GA] = double(NL_ZSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GB] = double(NL_ZSCALE(nl)) /
					  double(NL_XSCALE(nl))
	    do i = 1, NL_NPAR(nl) / 3 {
	        Memd[NL_PARSCALE(nl)+NL_GAMPL(i)] = double(NL_ZSCALE(nl))
	        Memd[NL_PARSCALE(nl)+NL_GCENT(i)] = double(NL_XSCALE(nl))
	        Memd[NL_PARSCALE(nl)+NL_GFWHM(i)] = double(NL_XSCALE(nl))
	    }

	case CGAUSS:
	    Memd[NL_PARSCALE(nl)+NL_GA]       = double(NL_ZSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GB]       = double(NL_ZSCALE(nl)) /
					        double(NL_XSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GAMPL(1)] = double(NL_ZSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GCENT(1)] = double(NL_XSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_GFWHM(1)] = double(NL_XSCALE(nl))
	    if (NL_NPAR(nl) > 5) {
	        do i = 2, NL_NPAR(nl) / 3 {
	            Memd[NL_PARSCALE(nl)+NL_GAMPL(i)]   = 1.d0
	            Memd[NL_PARSCALE(nl)+NL_GCENT(i)] = double(NL_XSCALE(nl))
	            Memd[NL_PARSCALE(nl)+NL_GFWHM(i)] = double(NL_XSCALE(nl))
	        }
	    }

	case USER:
	    do i = 0, NL_NPAR(nl)-1
	        Memd[NL_PARSCALE(nl)+i]   = 1.0d0	# No scaling

	case TWODGAUSS:
	    Memd[NL_PARSCALE(nl)+NL_G2A]    = double(NL_ZSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_G2AMPL] = double(NL_ZSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_G2XC]   = double(NL_XSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_G2YC]   = double(NL_YSCALE(nl))
	    Memd[NL_PARSCALE(nl)+NL_G2ELL]  = double(1.)
	    Memd[NL_PARSCALE(nl)+NL_G2TETA] = double(1.)
	    Memd[NL_PARSCALE(nl)+NL_G2FWHM] = double(NL_XSCALE(nl))

	default:
	    call error (0, "Unknown function.")
	}

	do i = 0, NL_NPAR(nl)-1
	    Memr[NL_SPARAMS(nl)+i] = Memr[NL_PARAMS(nl)+i] * 
				     Memd[NL_PARSCALE(nl)+i]
end

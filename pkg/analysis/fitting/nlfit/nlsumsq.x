include "nlfit.h"

define	PA		Memr[NL_SPARAMS(nl)+NL_G2TETA]
define	EL		Memr[NL_SPARAMS(nl)+NL_G2ELL]
define	PI		3.141592654
define	PI2		(PI / 2.)
define	DPI		(2 * PI)

# NL_SUMSQ --  Computes unweighted sum of squares and chi-square (both unscaled)
# between data and proposed fitting function. The returned function value is 
# the chi-square, which is the function to be minimized. This function is 
# called by the non-linear function minimization routine nl_amoeba. 
# Array NL_PFLAGS(nl) must be set to indicate which coefficients are being 
# varied by amoeba. Array NL_REJFLAG(nl) must be set to indicate which data 
# points are being rejected from the present fit. Chi-square computation only
# makes sense if weights are the data point's inverse errors.

real procedure nl_sumsq (nl, coeff)

pointer	nl		# i: Curve descriptor.
real	coeff[ARB]	# i: Array with coefficients which are being varyed
			#    by amoeba code.

#--
pointer	coeff1, sp
int	i, j, k
real	data, func, weig, resid
double	rms, chisq

int	nl_stati()

errchk	nl_eval

begin
	# Allocate area for coefficient array.
	call smark (sp)
	call salloc (coeff1, NL_NPAR(nl), TY_REAL)

	# Put in coefficient array both coefficients which are being varied
	# by amoeba code, as well as those which are fixed.
	j = 0
	do i = 0, NL_NPAR(nl)-1 {
	    if (Memb[NL_PFLAGS(nl)+i]) {
	        j = j + 1
	        # The following code is a patch to take care of off-limits
	        # values for ellipticity and position angle in 
	        # 2-D Gaussian fitting. It should be substituted
	        # by a more general approach to the problem of
	        # constrained fitting.
	        if (nl_stati (nl, "fitfunc") == TWODGAUSS) {
	            if (i == NL_G2ELL) {
	                if (coeff[j] < 0.)
	                    coeff[j] = -coeff[j]
	                if (coeff[j] > 1.)
	                    coeff[j] = 1 - coeff[j]
	            } else if (i == NL_G2TETA) {
	                while (coeff[j] > DPI)
	                    coeff[j] = coeff[j] - DPI
	                if (coeff[j] > PI)
	                    coeff[j] = coeff[j] - PI
	                while (coeff[j] < -DPI)
	                    coeff[j] = coeff[j] + DPI
	                if (coeff[j] < -PI)
	                    coeff[j] = coeff[j] + PI
	            }
	        }
		Memr[coeff1+i] = coeff[j]
	    }
	    else
		Memr[coeff1+i] = Memr[NL_SPARAMS(nl)+i]
	}

	# Obtain function array.
	call nl_eval (nl, Memr[NL_SX(nl)], Memr[NL_SY(nl)], 
                      Memr[NL_AUX(nl)], Memr[coeff1], NL_NPTS(nl))

	# Calculate residuals and sums, rejecting flagged points.
	rms   = 0.d0
	chisq = 0.d0
	k     = 0
	do i = 0, NL_NPTS(nl)-1 {
	    if ( !( Memb[NL_REJFLAG(nl)+i])) {
	        weig = Memr[NL_SW(nl)+i]  * NL_WSCALE(nl)
	        if (weig != 0.) {
	            data  = Memr[NL_SZ(nl)+i]  / NL_ZSCALE(nl)
	            func  = Memr[NL_AUX(nl)+i] / NL_ZSCALE(nl)
	            resid = data - func
		    rms   = rms + double(resid) * double(resid)
	            resid = resid * weig
	            chisq = chisq + double(resid) * double(resid)
	            k     = k + 1
	        }
	    }
	}
	call sfree (sp)

	# drop the factor -2 which made data set with only ncoeff+2 data points
	# crash. 6/28/94 JC Hsu
	#return ((rms / real (k - j - 2)) / (ssw / real(k)))

	if (k > 1) {
	    NL_RMS(nl) = (real(rms) / real(k-1))
	    if (NL_RMS(nl) > 0.)
	        NL_RMS(nl) = sqrt(NL_RMS(nl))
	    else
	        NL_RMS(nl) = 0.0          # Zero is used instead of INDEF
	} else                            # to avoid crash in replication
	    NL_RMS(nl) = 0.0              # sampling.

	if (k > j)
	    return (real(chisq) / real (k - j))
	else
	    call error (0, "Not enough data points")
end

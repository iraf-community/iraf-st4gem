include <error.h>
include "nlfit.h"

# NL_FIT -- Main non-linear function fitting procedure.
# This module fits a non-linear function to a set of data points. 
# The data points are defined by three related vectors: 2 independant variable 
# and one dependent variable vector. One-dimensional functions need only the
# first and third vectors to be defined by the caller. Also, a weigth array 
# must be supplied. This weight array must be built from the inverses of data
# error bars for the chi-square to make sense. The code assumes that work 
# areas were already been created and initialized by the caller, through the 
# use of procedure nl_init.

procedure nl_fit (nl, xfit, yfit, zfit, wt)

pointer nl			# i:  Curve descriptor.
real	xfit[ARB]		# i:  Independant x variable array.
real	yfit[ARB]		# i:  Independant y variable array.
real	zfit[ARB]		# i:  Dependant variable array.
real	wt[ARB]			# io: Weight array. This might be modified
				#     to suit chi-square computation.
#--
int	i, j
pointer	sp, sx, sx2, par
int	nstart, errc
double	zrms			# Raw RMS Z value.
double	hold
real	trms, noise
real	srms, schi
real	ccdn, rn2
real	maxwt, z, w
long	seed

real	nl_ngauss()
int	nl_minimize()

errchk	nl_minimize, nl_eval, nl_xscale

begin
	# Compute scaling factors for input variables.
	call nl_xscale (nl, xfit, yfit, NL_NPTS(nl))

	# Scale dependant variable and weights. Do not scale
        # if user function or Levenberg-Marquardt method.
	if (NL_FITFUNC(nl) == USER       ||
            NL_METHOD(nl)  == MARQUARDT) {

	    NL_ZSCALE(nl) = 1.0
	    NL_WSCALE(nl) = 1.0

	} else {

	    zrms = 0.
	    do i = 1, NL_NPTS(nl)
	        zrms = zrms + double(zfit[i])**2
	    zrms = sqrt(zrms/2.)
	    i = -log10(zrms)
	    NL_ZSCALE(nl) = 10.0d0 ** i
	    maxwt = 0.
	    do i = 1, NL_NPTS(nl)
	        maxwt = max(maxwt, abs(wt[i]))
	    if (maxwt <= 0.)
	        call error (0, "All weigths zero.")
	    else
	        NL_WSCALE(nl) = maxwt
	}

	# Realloc memory for fit.
	call realloc (NL_SX(nl), NL_NPTS(nl), TY_REAL)
	call realloc (NL_SY(nl), NL_NPTS(nl), TY_REAL)
	call realloc (NL_SZ(nl), NL_NPTS(nl), TY_REAL)
	call realloc (NL_SW(nl), NL_NPTS(nl), TY_REAL)

	# Generate scaled variables.
	call amulkr (xfit, NL_XSCALE(nl), Memr[NL_SX(nl)], NL_NPTS(nl))
	call amulkr (yfit, NL_YSCALE(nl), Memr[NL_SY(nl)], NL_NPTS(nl))
	call amulkr (zfit, NL_ZSCALE(nl), Memr[NL_SZ(nl)], NL_NPTS(nl))
	call aabsr  (wt, Memr[NL_SW(nl)], NL_NPTS(nl))
	call adivkr (Memr[NL_SW(nl)], NL_WSCALE(nl), Memr[NL_SW(nl)],
	             NL_NPTS(nl))

	# Now scale input parameters (user's guess).
	call nl_scale (nl)

	# Do it.
	nstart = 0
	errc = ERR
	while ((errc == ERR)              && 
               (nstart <= NL_RESTART(nl))) {
	    errc = nl_minimize (nl)
	    if ((errc == ERR) && (nstart < NL_RESTART(nl))) {
	        call eprintf ("  Starting again... ")
	        call flush (STDERR)
	    }
	    nstart = nstart + 1
	}

	# Get back to unscaled parameter values.
	do i = 0, NL_NPAR(nl)-1 {
	    Memr[NL_PARAMS(nl)+i]  = Memr[NL_SPARAMS(nl)+i] / 
				     Memd[NL_PARSCALE(nl)+i]
	    Memd[NL_PERRORS(nl)+i] = Memd[NL_PERRORS(nl)+i] / 
				     Memd[NL_PARSCALE(nl)+i]
	}

	# Compute coefficient errors by bootstrap resampling.
	if (NL_ERR(nl) && (NL_RMS(nl) > 0.)) {
	    seed = 1
	    call smark (sp)
	    call salloc (par, NL_NPAR(nl), TY_REAL)
	    call salloc (sx,  NL_NPAR(nl), TY_DOUBLE)
	    call salloc (sx2, NL_NPAR(nl), TY_DOUBLE)

	    # Initialize sums.
	    call amovr  (Memr[NL_SPARAMS(nl)], Memr[par], NL_NPAR(nl))
	    do i = 0, NL_NPAR(nl)-1 {
		hold = Memr[NL_SPARAMS(nl)+i] / Memd[NL_PARSCALE(nl)+i]
		Memd[sx+i]  = hold
		Memd[sx2+i] = hold * hold
	    }

	    # Store for posterior recovery at the end.
	    srms = NL_RMS(nl)
	    schi = NL_CHISQ(nl)

	    # Replicas are built in scaled space.
	    trms = NL_RMS(nl) * NL_ZSCALE(nl)

	    rn2 = NL_RNOISE(nl) ** 2  # ccd readout variance

	    do i = 1, NL_REPLI(nl) {

		# Generate replication = original fit curve + noise.
	        # Choose noise type from NL_ERRTYP(nl).
		call nl_eval (nl, Memr[NL_SX(nl)], Memr[NL_SY(nl)], 
			      Memr[NL_SZ(nl)], Memr[par], NL_NPTS(nl))

	        do j = 0, NL_NPTS(nl)-1 {
	            z = Memr[NL_SZ(nl)+j] 
	            w = Memr[NL_SW(nl)+j]

	            switch (NL_ERRTYP(nl)) {
	            case UNIFORM:
	                # Noise amplitude is own fit's rms.
	                noise = nl_ngauss (trms, seed)
	            case BARS:
	                # Noise amplitude is given by each datum error bar.
	                if (w <= 0.)
	                    noise = nl_ngauss (trms, seed)
	                else
	                    noise = nl_ngauss (NL_ZSCALE(nl) / 
	                                       NL_WSCALE(nl) / w, seed)
	            case CCD:
	                if (NL_EPADU(nl) > 0.0) {
	                    if (z <= (NL_ZSCALE(nl)/NL_EPADU(nl)))
	                        # one or less count.
	                        noise = nl_ngauss (NL_RNOISE(nl) * 
	                                           NL_ZSCALE(nl), seed)
	                    else
	                        # z must be translated to unscaled counts and
	                        # deviate scaled back. Use gaussian deviate
	                        # because poisson produces too large errors
	                        # (always ~20% of coefficient, no matter what 
	                        # the signal-to-noise ratio is). This was 
                                # veryfied in extensive tests with artificial 
	                        # data (IB 06/03/96).
	                        ccdn = z / NL_ZSCALE(nl) * NL_EPADU(nl) + rn2
	                        noise = nl_ngauss (sqrt(ccdn), seed) *
	                                NL_ZSCALE(nl)
	                } else 
	                    noise = nl_ngauss (NL_RNOISE(nl) * NL_ZSCALE(nl), 
	                                       seed)
	            }
		    Memr[NL_SZ(nl)+j] = z + noise
	        }

		# Now fit replicated data.
	        nstart = 0
	        errc = ERR
	        while ((errc == ERR)              && 
                       (nstart <= NL_RESTART(nl))) {
	            errc = nl_minimize (nl)
	            if ((errc == ERR) && (nstart < NL_RESTART(nl))) {
	                call eprintf ("  Starting again... ")
	                call flush (STDERR)
	            }
	            nstart = nstart + 1
	        }

		# Update sums.
		do j = 0, NL_NPAR(nl)-1 {
		    hold = Memr[NL_SPARAMS(nl)+j] / Memd[NL_PARSCALE(nl)+j]
		    Memd[sx+j]  = Memd[sx+j]  + hold
		    Memd[sx2+j] = Memd[sx2+j] + hold * hold
		}
	    }

	    # Calculate coeficient dispersions.
	    do j = 0, NL_NPAR(nl)-1 {
		Memd[NL_PERRORS(nl)+j] = Memd[sx2+j] - Memd[sx+j] ** 2 /
                                                       real(NL_REPLI(nl)+1)
		if (Memd[NL_PERRORS(nl)+j] > 0.)
		    Memd[NL_PERRORS(nl)+j] = sqrt(Memd[NL_PERRORS(nl)+j] /
					     real (NL_REPLI(nl)))
		else
		    Memd[NL_PERRORS(nl)+j] = 0.d0
	        if (!( Memb[NL_PFLAGS(nl)+j]))    # Corrects possible
		    Memd[NL_PERRORS(nl)+j] = 0.d0  # rounding errors.
	    }

	    call sfree (sp)

	    # Restore values from original fit.`
	    NL_RMS(nl)   = srms
	    NL_CHISQ(nl) = schi
	}
end

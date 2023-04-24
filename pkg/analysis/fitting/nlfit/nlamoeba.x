include <mach.h>
include "nlfit.h"

define	FTOL	5.0e-6	# Convergency criterion applied to rtol.

define	P	Memr[p    + np*$1 + $2]		# Array elements
define	PR	Memr[pr   + $1]
define	PRR	Memr[prr  + $1]
define	PBAR	Memr[pbar + $1]
define	Y	Memr[y    + $1]

# NL_AMOEBA -- Downhill simplex minimization method.
#
# This module implements the AMOEBA algorithm from Numerical Recipes, 
# which minimizes a function by the downhill simplex method. It calls a 
# function, nl_sumsq, which returns the chi-squared of residuals between 
# the data and the found solution. 
#
# On input, array NL_SPARAMS(nl) contains an initial guess, and on return 
# the fitted coefficients. On return, NL_CHISQ(nl) contains the chi-squared, 
# and NL_RMS(nl) contains the r.m.s residual of the fit. 
#
# On input, array NL_PFLAGS(nl) must be set in order to flag (true) 
# coefficients which are allowed to vary during the fit. 
#
# In case of no convergency, ERR is returned and a warning message is 
# issued, otherwise OK is returned.

int procedure nl_amoeba (nl)

pointer nl		# i: Curve descriptor.

#--

int	np		# Number of coefficients allowed to vary during fit.
int	ilo		# Indices of highest, next-highest and lowest
int	ihi		# simplex vertices.
int	inhi
pointer	sp
pointer	p		# Storage for coordinates of simplex vertices.
pointer	y		# Storage for function values on simplex vertices.
pointer	pr, prr, pbar	# Auxiliary storage areas.
pointer	dev		# Deviations to build initial simplex.
real	ypr, yprr
real	rtol		# Relative function variation in current iteration.
pointer	lrtol		# Circular buffer for storing previous values of rtol.
int	irtol
bool	flag, store
int	iter		# Number of iterations.
int 	ret
int	i, ii, j, jj

real	nl_sumsq()	# Evaluates sum of squares.

errchk	nl_sumsq, salloc

begin
	# Set up number of variable coefficients.
	np = 0
	do i = 0, NL_NPAR(nl) - 1 {
	    if (Memb[NL_PFLAGS(nl)+i])
	        np = np + 1
	}
	if (np == 0)
	    call error (0, "Less than 1 parameter for current fit.")

	# Allocate storage for simplex-related quantities.
	call smark (sp)
	call salloc (p,    (np+1) * np, TY_REAL)
	call salloc (y,    np+1, TY_REAL)
	call salloc (pr,   np,   TY_REAL)
	call salloc (prr,  np,   TY_REAL)
	call salloc (pbar, np,   TY_REAL)

	call salloc (lrtol, np+1, TY_REAL)
	call amovkr (MAX_REAL, Memr[lrtol], np+1)

	# Create initial simplex. Uses only coefficients which are to
	# be varied. Deviations from supplied coefficients depend on 
	# functional form:
	# Gaussian types:          ampl + 5 %
	#                          cent + 5 % FWHM
	#                          FWHM + 5 %
	# other functions: coeff + 10 %

	call salloc (dev, NL_NPAR(nl), TY_REAL)

	if ((NL_FITFUNC(nl) == GAUSSIANS) || (NL_FITFUNC(nl) == CGAUSS)) {
	    Memr[dev+NL_GA] = 0.1 * Memr[NL_SPARAMS(nl) + NL_GA]
	    Memr[dev+NL_GB] = 0.1 * Memr[NL_SPARAMS(nl) + NL_GB]
	    do i = 1, ((NL_NPAR(nl)-2) / 3) {
	        Memr[dev+NL_GAMPL(i)] = 0.05*Memr[NL_SPARAMS(nl) + NL_GAMPL(i)]
	        Memr[dev+NL_GCENT(i)] = 0.05*Memr[NL_SPARAMS(nl) + NL_GFWHM(i)]
	        Memr[dev+NL_GFWHM(i)] = 0.05*Memr[NL_SPARAMS(nl) + NL_GFWHM(i)]
	    }
	} else if (NL_FITFUNC(nl) == TWODGAUSS) {
	    do i = 0, NL_NPAR(nl) - 1 
	        Memr[dev+i] = Memr[NL_SPARAMS(nl)+i] * 0.1  # 10% deviations.
	    Memr[dev+NL_G2AMPL] = 0.05*Memr[NL_SPARAMS(nl) + NL_G2AMPL]
	    Memr[dev+NL_G2XC]   = 0.05*Memr[NL_SPARAMS(nl) + NL_G2FWHM]
	    Memr[dev+NL_G2YC]   = 0.05*Memr[NL_SPARAMS(nl) + NL_G2FWHM]
	    Memr[dev+NL_G2FWHM] = 0.05*Memr[NL_SPARAMS(nl) + NL_G2FWHM]
	} else {
	    do i = 0, NL_NPAR(nl) - 1 
	        Memr[dev+i] = Memr[NL_SPARAMS(nl)+i] * 0.1  # 10% deviations.
	}

	ii = 0
	do i = 0, NL_NPAR(nl) - 1 {
	    if (Memb[NL_PFLAGS(nl)+i]) {
	        jj = 0
	        do j = 0, NL_NPAR(nl) - 1 {
	            if (Memb[NL_PFLAGS(nl)+j])  {
	                P(ii,jj) = Memr[NL_SPARAMS(nl)+j]
	                jj = jj + 1
	            }
	        }
	        ii = ii + 1
	    }
	}
	jj = 0
	do j = 0, NL_NPAR(nl) - 1 {
	    if (Memb[NL_PFLAGS(nl)+j])  {
	        P(np,jj) = Memr[NL_SPARAMS(nl)+j]
	        jj = jj + 1
	    }
	}

	# Apply deviations to build simplex vertices.
	jj = 0
	do j = 0, NL_NPAR(nl) - 1 {
	    if (Memb[NL_PFLAGS(nl)+j])  {
	        P(jj,jj) = P(jj,jj) + Memr[dev+j]
	        jj = jj + 1
	    }
	}

	# Get initial chi-square values.
	do i = 0, np {
	    do j = 0, np-1
	        PR(j) = P(i,j)
	    Y(i) = nl_sumsq(nl, Memr[pr])
	}


# .................... Start of AMOEBA code ..............................

	store = false
	irtol = 0
	iter  = 0
	repeat {
	    ilo = 0
	    if (Y(0) > Y(1)) {
	        ihi = 0
	        inhi = 1
	    } else {
	        ihi = 1
	        inhi = 0
	    }
	    do i = 0, np {
 	        if (Y(i) <= Y(ilo)) 
	            ilo = i
	        if (Y(i) > Y(ihi)) {
	            inhi = ihi
	            ihi = i
	        } else if (Y(i) > Y(inhi)) {
	            if (i != ihi)
	                inhi = i
	        }
	    }

	    rtol = 2. * abs(Y(ihi) - Y(ilo)) / (abs(Y(ihi)) + abs(Y(ilo)))

	    if (NL_VERB(nl)) {
	        call eprintf ("iter= %d  rtol= %g\n")
	            call pargi (iter)
	            call pargr (rtol)
	            call flush (STDERR)
	    }

	    # Converged !
	    if (rtol < FTOL) {
	        ret = OK
	        break
	    }
	    # Did not converge, check last np+1 remembered iterations.
	    store = ! store
	    if (store) {
	        Memr[lrtol+irtol] = rtol
	        irtol = irtol + 1
	        if (irtol > np)
	            irtol = 0
	    }
	    flag = true
	    do i = 0, np -1 {
	        do j = i+1, np {
	            if (Memr[lrtol+i] != Memr[lrtol+j])
                        flag = false
	        }
	    }
	    # Issue msg and stop.
	    if ((iter == NL_MAXIT(nl)) || (flag)) {
	        call eprintf("\nStop after %d iterations. Tolerance = %g  ")
	            call pargi (iter)
	            call pargr (rtol)
	        call flush (STDERR)
	        ret = ERR
	        break
	    }

	    iter = iter + 1
	    do j = 0, np-1
	        PBAR(j) = 0.
	    do i = 0, np {
	        if (i != ihi) {
	            do j = 0, np-1
	                PBAR(j) = PBAR(j) + P(i,j)
	        }
	    }
	    do j = 0, np-1 {
	        PBAR(j) = PBAR(j) / np
	        PR(j) = (1. + NL_ALPHA(nl)) * PBAR(j) - NL_ALPHA(nl) * P(ihi,j)
	    }
	    ypr = nl_sumsq(nl, Memr[pr])
	    if (ypr < Y(ilo)) {
	        do j = 0, np-1
	            PRR(j) = NL_GAMMA(nl) * PR(j) + (1. - NL_GAMMA(nl)) * 
	                     PBAR(j)
	        yprr = nl_sumsq(nl, Memr[prr])
	        if (yprr < Y(ilo)) {
	            do j = 0, np-1
	                P(ihi,j) = PRR(j)
	            Y(ihi) = yprr
	        } else {
	            do j = 0, np-1
	                P(ihi,j) = PR(j)
	            Y(ihi) = ypr
	        }
	    } else if (ypr >= Y(inhi)) {
	        if (ypr < Y(ihi)) {
	            do j = 0, np-1
	                P(ihi,j) = PR(j)
	            Y(ihi) = ypr
	        }
	        do j = 0, np-1
	            PRR(j) = NL_BETA(nl) * P(ihi,j) + (1.-NL_BETA(nl)) * PBAR(j)
	        yprr = nl_sumsq(nl, Memr[prr])
	        if (yprr < Y(ihi)) {
	            do j = 0, np-1
	                P(ihi,j) = PRR(j)
	            Y(ihi) = yprr
	        } else {
	            do i = 0, np {
	                if (i != ilo) {
	                    do j = 0, np-1 {
	                        PRR(j) = 0.5 * (P(i,j) + P(ilo,j))
	                        P(i,j) = PR(j)
	                    }
	                    Y(i) = nl_sumsq(nl, Memr[prr])
	                }
	            }
	        }
	    } else {
	        do j = 0, np-1
	            P(ihi,j) = PR(j)
	        Y(ihi) = ypr
	    }
	}

# ......................  End of AMOEBA code ...........................

#     Solution is in arrays p and y; transfer it to output array NL_SPARAMS(nl)
#     and NL_CHISQ(nl) parameter.

	jj = 0
	do j = 0, NL_NPAR(nl) - 1 {
	    if (Memb[NL_PFLAGS(nl)+j]) {
	        Memr[NL_SPARAMS(nl)+j] = P(1,jj)
	        jj = jj + 1
	    }
	}
	NL_CHISQ(nl) = Y(0)

	call sfree (sp)
	return (ret)
end

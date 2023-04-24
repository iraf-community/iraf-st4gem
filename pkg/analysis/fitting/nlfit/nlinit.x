include "nlfit.h"

# NL_INIT --  Initialize the curve structure and allocate memory 
# for the curve-related parameters. Routine assumes that structure 
# pointer has the value NULL when called for the first time.

procedure nl_init (nl, fitfunc, coeff, coefflags, npar, npix)

pointer	nl		# io: curve structure pointer
int	fitfunc		# i:  type of function to be fitted
real	coeff[ARB]	# i:  initial guess for coefficients
bool	coefflags[ARB]	# i:  coefficient flags
int	npar		# i:  number of coefficients
int	npix		# i:  number of data points

#--
int	i

errchk	realloc, malloc

begin
#	if (npar < 1)
#	    call error (0, "Illegal number of parameters.")
	if (npix < npar+1)
	    call error (0, "Insufficient number of data points.")
#	do i = 1, npar {
#	    if (IS_INDEFR (coeff[i]))
#	        call error (0, "There are undefined coefficients.")
#	}

	# Allocate space for the curve structure
	if (nl == NULL) {

	    call malloc (nl, LEN_NLSTRUCT, TY_STRUCT)

	    # Enforce allocation of brand new areas by realloc.
	    NL_PARAMS(nl)	= NULL
	    NL_SPARAMS(nl)	= NULL
	    NL_PARSCALE(nl)	= NULL
	    NL_PFLAGS(nl)	= NULL
	    NL_PERRORS(nl)	= NULL
	    NL_SX(nl)		= NULL
	    NL_SY(nl)		= NULL
	    NL_SZ(nl)		= NULL
	    NL_SW(nl)		= NULL
	    NL_AUX(nl)		= NULL
	    NL_REJFLAG(nl)	= NULL
	    NL_USERFUNC(nl)	= NULL
	    NL_UCODE(nl)	= NULL

	    # Set to default values.
	    NL_YSCALE(nl)   = 1.
	    NL_XSCALE(nl)   = 1.
	    NL_ZSCALE(nl)   = 1.
	    NL_WSCALE(nl)   = 1.
	    NL_ALPHA(nl)    = 1.0
	    NL_BETA(nl)     = 0.5
	    NL_GAMMA(nl)    = 2.0
	    NL_METHOD(nl)   = AMOEBA
	    NL_SIGMA(nl)    = INDEFR
	    NL_RMS(nl)	    = INDEFR
	    NL_CHISQ(nl)    = INDEFR
	    NL_UNIT(nl)     = AUTO
	    NL_ERR(nl)      = false
	    NL_ERRTYP(nl)   = UNIFORM
	    NL_MAXIT(nl)    = MAXITER
	    NL_RESTART(nl)  = DEF_RESTART
	    NL_REPLI(nl)    = DEF_REPLIC
	    NL_VERB(nl)     = false
	}

	# Set parameters for present fit.
	NL_FITFUNC(nl) = fitfunc
	IF (fitfunc >= NLTWODSEP)
	    NL_DIM(nl) = 2
	else
	    NL_DIM(nl) = 1
	NL_NPTS(nl)    = npix
	NL_NPAR(nl)    = npar

	# Allocate areas for curve-associated vectors and initialize.
	call realloc (NL_PARAMS(nl),   NL_NPAR(nl), TY_REAL)
	call realloc (NL_SPARAMS(nl),  NL_NPAR(nl), TY_REAL)
	call realloc (NL_PARSCALE(nl), NL_NPAR(nl), TY_DOUBLE)
	call realloc (NL_PFLAGS(nl),   NL_NPAR(nl), TY_BOOL)
	call realloc (NL_PERRORS(nl),  NL_NPAR(nl), TY_DOUBLE)

	call amovr  (coeff, Memr[NL_PARAMS(nl)],   npar)
	call amovr  (coeff, Memr[NL_SPARAMS(nl)],  npar)
	call amovi  (coefflags, Memb[NL_PFLAGS(nl)], npar)
	call amovkd (1.d0,  Memd[NL_PARSCALE(nl)], npar)
	call amovkd (0.d0,  Memd[NL_PERRORS(nl)],  npar)

	# Allocate areas for data points.
	call realloc (NL_SX(nl),      npix, TY_REAL)
	call realloc (NL_SY(nl),      npix, TY_REAL)
	call realloc (NL_SZ(nl),      npix, TY_REAL)
	call realloc (NL_SW(nl),      npix, TY_REAL)
	call realloc (NL_AUX(nl),     npix, TY_REAL)
	call realloc (NL_REJFLAG(nl), npix, TY_BOOL)

	# No rejected points yet.
	NL_NREJECT(nl) = 0
	do i = 0, npix-1
	    Memb[NL_REJFLAG(nl)+i] = false
end


# NL_FREE -- Free memory used by curve structure and data areas.

procedure nl_free (nl)

pointer	nl	# i: the curve descriptor

#--
errchk	mfree

begin
	if (nl == NULL)
	    return

	if (NL_SX(nl) != NULL) {
		call mfree (NL_SX(nl), TY_REAL)
		call mfree (NL_SY(nl), TY_REAL)
		call mfree (NL_SZ(nl), TY_REAL)
		call mfree (NL_SW(nl), TY_REAL)
		call mfree (NL_AUX(nl), TY_REAL)
		call mfree (NL_REJFLAG(nl), TY_BOOL)
	}

	if (NL_PARAMS(nl) != NULL) {
		call mfree (NL_PARAMS(nl), TY_REAL)
		call mfree (NL_SPARAMS(nl), TY_REAL)
		call mfree (NL_PARSCALE(nl), TY_DOUBLE)
		call mfree (NL_PFLAGS(nl), TY_BOOL)
		call mfree (NL_PERRORS(nl), TY_DOUBLE)
	}

	if (NL_USERFUNC(nl) != NULL) {
	        call mfree (NL_USERFUNC(nl), TY_CHAR)
		call nl_freec (nl)
	}

	call mfree (nl, TY_STRUCT)
end

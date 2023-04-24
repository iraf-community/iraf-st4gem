include "nlfit.h"

# This file contains procedures for seting and retrieving values
# from the curve structure.

# NL_STATI -- Return integer parameters from the nlfit package.

int procedure nl_stati (nl, param)

pointer	nl			# i: NLFIT pointer
char	param[ARB]		# i: Parameter to be gotten

#--
bool	streq()

begin
	if (streq (param, "fitfunc"))
	    return (NL_FITFUNC(nl))
	else if (streq (param, "algorithm"))
	    return (NL_METHOD(nl))
	else if (streq (param, "npar"))
	    return (NL_NPAR(nl))
	else if (streq (param, "npts"))
	    return (NL_NPTS(nl))
	else if (streq (param, "nreject"))
	    return (NL_NREJECT(nl))
	else if (streq (param, "units"))
	    return (NL_UNIT(nl))
	else if (streq (param, "maxit"))
	    return (NL_MAXIT(nl))
	else if (streq (param, "replic"))
	    return (NL_REPLI(nl))
	else if (streq (param, "restart"))
	    return (NL_RESTART(nl))
	else if (streq (param, "errtyp"))
	    return (NL_ERRTYP(nl))
	else if (streq (param, "dimension"))
	    return (NL_DIM(nl))
	else
	    call error (0, "Unknown integer parameter.")
end


# NL_STATB -- Return boolean parameters from the nlfit package.

bool procedure nl_statb (nl, param)

pointer	nl			# i: NLFIT pointer
char	param[ARB]		# i: Parameter to be gotten

#--
bool	streq()

begin
	if (streq (param, "errors"))
	    return (NL_ERR(nl))
	if (streq (param, "verbose"))
	    return (NL_VERB(nl))
	else
	    call error (0, "Unknown boolean parameter")
end


# NL_STATR -- Return real parameters from the nlfit package.

real procedure nl_statr (nl, param)

pointer	nl			# i: NLFIT pointer
char	param[ARB]		# i: Parameter to be gotten

#--
bool	streq()

begin
	if (streq(param, "rms"))
	    return (NL_RMS(nl))
	else if (streq(param, "chisq"))
	    return (NL_CHISQ(nl))
	else if (streq(param, "sigma"))
	    return (NL_SIGMA(nl))
	else if (streq(param, "epadu"))
	    return (NL_EPADU(nl))
	else if (streq(param, "readnoise"))
	    return (NL_RNOISE(nl))
	else if (streq(param, "xscale"))
	    return (NL_XSCALE(nl))
	else if (streq(param, "yscale"))
	    return (NL_YSCALE(nl))
	else if (streq(param, "zscale"))
	    return (NL_ZSCALE(nl))
	else if (streq(param, "wscale"))
	    return (NL_WSCALE(nl))
	else if (streq(param, "alpha"))
	    return (NL_ALPHA(nl))
	else if (streq(param, "beta"))
	    return (NL_BETA(nl))
	else if (streq(param, "gamma"))
	    return (NL_GAMMA(nl))
	else 
	    call error (0, "Unknown integer parameter.")
end


# NL_PUTI -- Put integer valued parameters.

procedure nl_puti (nl, param, ival)

pointer	nl			# i: NLFIT pointer
char	param[ARB]		# i: Parameter to be put
int	ival			# i: Integer value

#--
bool	streq()

begin
	if (streq (param, "fitfunc"))
	    NL_FITFUNC(nl) = ival
	else if (streq (param, "algorithm"))
	    NL_METHOD(nl) = ival
	else if (streq (param, "npar"))
	    NL_NPAR(nl) = ival
	else if (streq (param, "npts"))
	    NL_NPTS(nl) = ival
	else if (streq (param, "nreject"))
	    NL_NREJECT(nl) = ival
	else if (streq (param, "units"))
	    NL_UNIT(nl) = ival
	else if (streq (param, "errtyp"))
	    NL_ERRTYP(nl) = ival
	else if (streq (param, "maxit"))
	    NL_MAXIT(nl) = ival
	else if (streq (param, "replic"))
	    NL_REPLI(nl) = ival
	else if (streq (param, "restart"))
	    NL_RESTART(nl) = ival
	else if (streq (param, "dimension"))
	    NL_DIM(nl) = ival
	else
	    call error (0, "Unknown integer parameter.")
end


# NL_PUTB -- Put boolean valued parameters.

procedure nl_putb (nl, param, bval)

pointer	nl			# i: NLFIT pointer
char	param[ARB]		# i: Parameter to be put
bool	bval			# i: Boolean value

#--
bool	streq()

begin
	if (streq (param, "errors"))
	    NL_ERR(nl) = bval
	else if (streq (param, "verbose"))
	    NL_VERB(nl) = bval
	else
	    call error (0, "Unknown boolean parameter.")
end


# NL_PUTR -- Put real valued parameters.

procedure nl_putr (nl, param, rval)

pointer	nl			# i: NLFIT pointer
char	param[ARB]		# i: Parameter to be put
real	rval			# i: Real value

#--
bool	streq()

begin
	if (streq (param, "rms"))
	    NL_RMS(nl) = rval
	else if (streq (param, "chisq"))
	    NL_CHISQ(nl) = rval
	else if (streq (param, "sigma"))
	    NL_SIGMA(nl) = rval
	else if (streq (param, "epadu"))
	    NL_EPADU(nl) = rval
	else if (streq (param, "readnoise"))
	    NL_RNOISE(nl) = rval
	else if (streq (param, "xscale"))
	    NL_XSCALE(nl) = rval
	else if (streq (param, "yscale"))
	    NL_YSCALE(nl) = rval
	else if (streq (param, "zscale"))
	    NL_ZSCALE(nl) = rval
	else if (streq (param, "wscale"))
	    NL_WSCALE(nl) = rval
	else if (streq (param, "alpha"))
	    NL_ALPHA(nl) = rval
	else if (streq (param, "beta"))
	    NL_BETA(nl) = rval
	else if (streq (param, "gamma"))
	    NL_GAMMA(nl) = rval
	else
	    call error (0, "Unknown real parameter.")
end


# NL_GCOEFF -- Get the number and magnitude of the coefficients
# from the curve descriptor.

procedure nl_gcoeff (nl, coeff, npar)

pointer	nl		# i: curve descriptor
real	coeff[ARB]	# o: the coefficients of the fit
int	npar		# o: the number of coefficients

#--
begin
	npar = NL_NPAR(nl)
	call amovr (Memr[NL_PARAMS(nl)], coeff, npar)
end


# NL_GERRORS -- Get the number and magnitude of the coefficient errors 
# from the curve descriptor.

procedure nl_gerrors (nl, perror, npar)

pointer	nl		# i: curve descriptor
real	perror[ARB]	# o: errors
int	npar		# o: number of coefficients

#--
int	i

begin
	npar = NL_NPAR(nl)
	do i = 1, NL_NPAR(nl)
	    perror[i] = Memd[NL_PERRORS(nl)+i-1]
end


# NL_GFLAGS -- Get the coefficient flags.

procedure nl_gflags (nl, flags, npar)

pointer	nl		# i: curve descriptor
bool	flags[ARB]	# o: coefficient flags
int	npar		# o: number of coefficients

#--
begin
	npar = NL_NPAR(nl)
	call amovi (Memb[NL_PFLAGS(nl)], flags, npar)
end


# NL_SCOEFF -- Set the number and magnitude of the coefficients
# in the curve descriptor.

procedure nl_scoeff (nl, coeff, npar)

pointer	nl		# i: curve descriptor
real	coeff[ARB]	# i: the coefficients of the fit
int	npar		# i: the number of coefficients

#--
begin
	NL_NPAR(nl) = npar
	call amovr (coeff, Memr[NL_PARAMS(nl)], npar)
end


# NL_SERRORS -- Restore the error vector in the curve descriptor.

procedure nl_serrors (nl, errors, npar)

pointer	nl		# i: curve descriptor
real	errors[ARB]	# i: errors to be stored
int	npar		# i: the number of coefficients

#-
int	i

begin
	NL_NPAR(nl) = npar
	do i = 1, NL_NPAR(nl)
	    Memd[NL_PERRORS(nl)+i-1] = errors[i]
end


# NL_SFLAGS -- Set the coefficient flags.

procedure nl_sflags (nl, flags, npar)

pointer	nl		# i: curve descriptor
bool	flags		# i: flags array
int	npar		# i: number of coefficients

#--
begin
	call amovi (flags, Memb[NL_PFLAGS(nl)], npar)
end


# NL_IUSER  --  Initializes and compiles the user function string.

procedure nl_iuser (nl, function, size)

pointer	nl		# i: curve descriptor
char	function[ARB]	# i: string with function
int	size		# i: string size

#--

errchk  nl_compile

begin
	if (NL_FITFUNC(nl) == USER) {
	    call realloc (NL_USERFUNC(nl), size, TY_CHAR)
	    call strcpy (function, Memc[NL_USERFUNC(nl)], size)
	    call nl_compile (nl, function)
	}
end

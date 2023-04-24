include	<gset.h>
include	<error.h>
include	"icfit.h"

# List of colon commands.
define	CMDS "|show|sample|naverage|function|low_reject|high_reject\
	|niterate|grow|chisq|rms|vshow|coef|fix|var|resample|maxit|verbose\
        |ngauss|pcomp|xaxis|yaxis|mag0|ltype|alpha|beta|gamma|user\
	|replic|restart|method"

define	SHOW		1	# Show values of parameters
define	SAMPLE		2	# Set or show sample ranges
define	NAVERAGE	3	# Set or show sample averaging or medianing
define	FUNCTION	4	# Set or show function type
define	LOW_REJECT	5	# Set or show lower rejection factor
define	HIGH_REJECT	6	# Set or show upper rejection factor
define	NITERATE	7	# Set or show rejection iterations
define	GROW		8	# Set or show rejection growing radius
define	CHISQ		9	# Show chisq of fit
define	RMS		10	# Show rms of fit
define	VSHOW		11	# Show verbose information
define  COEF		12	# Set or show function coefficient
define	FIX		13	# Fix coefficient
define	VARY		14	# Vary coefficient
define	ERRORS		15	# Set or show coeff. error evaluation
define	MAXIT		16	# Set or show maximum number of iterations
define	VERBOSE		17	# Set or show iteration info status
define	NGAUSS		18	# Set or show number of Gaussians
define	PCOMP		19	# Set or show component plotting.
define	XAXIS		20	# Set or show type of x axis
define	YAXIS		21	# Set or show type of y axis
define	MAG0		22	# Set or show magnitude axis constant
define	LTYPE		23	# Set or show line type (continuous|boxes)
define	ALPHA		24	# Set or show alpha
define	BETA		25	# Set or show beta
define	GAMMA		26	# Set or show gamma
define	USF		27	# Set or show user function string
define	REPLIC		28	# Set or show # of resamplin replicas
define	RESTART		29	# Set or show # of restarts
define	METHOD		30	# Set or show current minimzation method


# ICG_COLON -- Processes cursor colon commands.  
# The common flags and newgraph signal changes in fitting parameters or the 
# need to redraw the graph.

procedure icg_colon (ic, cmdstr, newgraph, newaxis, gp, gt, nl, x, y, wts, npts)

pointer	ic				# ICFIT pointer
char	cmdstr[ARB]			# Command string
int	newgraph			# New graph?
int	newaxis				# New yaxis?
pointer	gp				# GIO pointer
pointer	gt				# GTOOLS pointer
pointer	nl				# NLFIT pointer
real	x[npts], y[npts], wts[npts]	# Data arrays
int	npts				# Number of data points
#--

char	cmd[SZ_LINE]
int	ncmd, ival, newf, newm, i
real	rval, rval2
bool	bval

int	nscan(), strdic(), icg_col2(), nl_stati()
real	nl_statr()
bool	nl_statb()

begin
	# Use formated scan to parse the command string.
	# The first word is the command and it may be minimum match
	# abbreviated with the list of commands.

	call sscan (cmdstr)
	call gargwrd (cmd, SZ_LINE)
	ncmd = strdic (cmd, cmd, SZ_LINE, CMDS)

	switch (ncmd) {
	case SHOW: # :show - Show the values of the fitting parameters.
	    call gargwrd (cmd, SZ_LINE)
	    if (nscan() == 1) {
		call gdeactivate (gp, AW_CLEAR)
		call ic_show (ic, "STDOUT", gt, nl)
		call greactivate (gp, AW_PAUSE)
	    } else {
		iferr (call ic_show (ic, cmd, gt, nl))
		    call erract (EA_WARN)
	    }

	case SAMPLE: # :sample - List or set the sample points.
	    call gargstr (cmd, SZ_LINE)
	    if (cmd[1] == EOS) {
	        call printf ("sample = %s\n")
		    call pargstr (Memc[IC_SAMPLE(ic)])
	    } else {
		call strcpy (cmd, Memc[IC_SAMPLE(ic)], SZ_LINE)
		IC_NEWX(ic) = YES
	    }

	case NAVERAGE: # :naverage - List or set the sample averging.
	    call gargi (ival)
	    if (nscan() == 1) {
		call printf ("naverage = %d\n")
		    call pargi (IC_NAVERAGE(ic))
	    } else {
		IC_NAVERAGE(ic) = ival
		IC_NEWX(ic) = YES
	    }

	case FUNCTION: # :function - List or set the fitting function.
	    call gargwrd (cmd, SZ_LINE)
	    if (nscan() == 1) {
		call printf ("function = %s\n")
		    call ic_gstr (ic, "function", cmd, SZ_LINE)
		    call pargstr (cmd)
	    } else {
	        newf = strdic (cmd, cmd, SZ_LINE, FUNCTIONS)
	        if (newf != 0) {
	            call icg_col4 (ic, nl, x, y, wts, newf, 0) 
		    IC_NEWFUNCTION(ic) = YES
	        } else 
	            call eprintf ("Non-recognizable name.\n")
	    }

	case USF: # :user - List or set the user fitting function.
	    if (nl_stati(nl, "fitfunc") == USER) {
	        call gargwrd (cmd, SZ_LINE)
	        if (nscan() == 1) {
		    call printf ("function = %s\n")
		        call pargstr (Memc[NL_USERFUNC(nl)])
	        } else
	            call nl_iuser (nl, cmd, SZ_LINE)
	    }

	case XAXIS: # :xaxis - List or set the x axis type.
	    call gargwrd (cmd, SZ_LINE)
	    if (nscan() == 1) {
		call printf ("xaxis = %s\n")
	        switch (IC_XAXIS(ic)) {
	        case IC_LOG:
	            call pargstr ("log")
	        default:
	            call pargstr ("linear")
	        }
	    } else {
	        i = strdic (cmd, cmd, SZ_LINE, AXISTYPE)
	        if (i != 0) {
	            IC_XAXIS(ic) = i - 1
	            newgraph = YES
	            newaxis  = YES
	        } else 
	            call eprintf ("Non-recognizable name.\n")
	    }

	case YAXIS: # :yaxis - List or set the y axis type.
	    call gargwrd (cmd, SZ_LINE)
	    if (nscan() == 1) {
		call printf ("yaxis = %s\n")
	        switch (IC_YAXIS(ic)) {
	        case IC_LOG:
	            call pargstr ("log")
	        case IC_MAG:
	            call pargstr ("mag")
	        default:
	            call pargstr ("linear")
	        }
	    } else {
	        i = strdic (cmd, cmd, SZ_LINE, AXISTYPE)
	        if (i != 0) {
	            IC_YAXIS(ic) = i - 1
	            newgraph = YES
	            newaxis  = YES
	        } else 
	            call eprintf ("Non-recognizable name.\n")
	    }

	case MAG0: # :mag0 - List or set the magnitude scale constant.
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("mag0 = %g\n")
		    call pargr (IC_MAG0(ic))
	    } else
		IC_MAG0(ic) = rval

	case LTYPE: # :ltype - select type of plot symbol
	    call gargwrd (cmd, SZ_LINE)
	    if (nscan() == 1) {
		call printf ("ltype = %s\n")
	        switch (IC_LTYPE(ic)) {
	        case IC_BOX:
	            call pargstr ("boxes")
	        default:
	            call pargstr ("continuous")
	        }
	    } else {
	        i = strdic (cmd, cmd, SZ_LINE, "|continuous|boxes|")
	        if (i != 0) {
	            IC_LTYPE(ic) = i - 1
	            newgraph = YES
	            newaxis  = YES
	        } else 
	            call eprintf ("Non-recognizable name.\n")
	    }

	case LOW_REJECT: # :low_reject - List or set lower rejection factor.
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("low_reject = %g\n")
		    call pargr (IC_LOW(ic))
	    } else
		IC_LOW(ic) = rval

	case HIGH_REJECT: # :high_reject - List or set high rejection factor.
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("high_reject = %g\n")
		    call pargr (IC_HIGH(ic))
	    } else
		IC_HIGH(ic) = rval

	case NITERATE: # :niterate - List or set the rejection iterations.
	    call gargi (ival)
	    if (nscan() == 1) {
		call printf ("niterate = %d\n")
		    call pargi (IC_NITERATE(ic))
	    } else {
		IC_NITERATE(ic) = ival
	    }

	case GROW: # :grow - List or set the rejection growing.
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("grow = %g\n")
		    call pargr (IC_GROW(ic))
	    } else
		IC_GROW(ic) = rval

	case CHISQ: # :chisq - List fit's chi-sq.
	    call gargr (rval)
	    rval = nl_statr (nl, "chisq")
	    call printf ("chi-sq = %g\n")
	        call pargr (rval)

	case RMS: # :rms - List fit's rms.
	    call gargr (rval)
	    rval = nl_statr (nl, "rms")
	    call printf ("rms = %g\n")
	        call pargr (rval)

	case VSHOW: # :vshow - Verbose list of the fitting parameters. 
	    call gargwrd (cmd, SZ_LINE)
	    if (nscan() == 1) {
		call gdeactivate (gp, AW_CLEAR)
		call ic_vshow (ic, "STDOUT", nl, x, y, wts, npts, gt)
		call greactivate (gp, AW_PAUSE)
	    } else {
		iferr {
		    call ic_vshow (ic, cmd, nl, x, y, wts, npts, gt)
		} then 
		    call erract (EA_WARN)
	    }

	case COEF: # :coef - list or set function coefficient.
	    call gargi (ival)
	    call gargr (rval)
	    if (nscan() == 1) {
	        return
	    } else if (nscan() == 2) {
		if (icg_col2 (nl, ival, rval, rval2, bval) == OK) {
		    call printf ("coef%d = %g (%g)  %s\n")
	            call pargi (ival)
		    call pargr (rval)
	            call pargr (rval2)
	            if (bval)
	                call pargstr ("var")
                    else
                        call pargstr ("fix")
	        }
	    } else {
		call icg_col1 (nl, ival, rval)
		IC_NEWFUNCTION(ic) = YES
	    }

	case FIX: # :fix - fix coefficient(s)
	    call gargi (ival)
	    if (nscan() == 1)
	        ival = 0
	    call icg_col3 (nl, ival, false)
	    IC_NEWFUNCTION(ic) = YES

	case VARY: # :vary - vary coefficient(s)
	    call gargi (ival)
	    if (nscan() == 1)
	        ival = 0
	    call icg_col3 (nl, ival, true)
	    IC_NEWFUNCTION(ic) = YES

	case ERRORS: # :errors - List or set the error evaluation control flag
	    call gargb (bval)
	    if (nscan() == 1) {
		call printf ("resample = %b\n")
		    call pargb (nl_statb (nl, "errors"))
	    } else
		call nl_putb (nl, "errors", bval)

	case ALPHA: # :alpha - Set or list the amoeba alpha parameter
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("alpha = %g\n")
		    call pargr (nl_statr (nl, "alpha"))
	    } else
		call nl_putr (nl, "alpha", rval)

	case BETA: # :beta - Set or list the amoeba beta parameter
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("beta = %g\n")
		    call pargr (nl_statr (nl, "beta"))
	    } else
		call nl_putr (nl, "beta", rval)

	case GAMMA: # :gamma - Set or list the amoeba gamma parameter
	    call gargr (rval)
	    if (nscan() == 1) {
		call printf ("gamma = %g\n")
		    call pargr (nl_statr (nl, "gamma"))
	    } else
		call nl_putr (nl, "gamma", rval)

	case MAXIT: # :maxit - Set or list the maximum no. of iterations
	    call gargi (ival)
	    if (nscan() == 1) {
		call printf ("maxit = %d\n")
		    call pargi (nl_stati (nl, "maxit"))
	    } else
		call nl_puti (nl, "maxit", ival)

	case REPLIC: # :replic - List or set # of replicas.
	    call gargi (ival)
	    if (nscan() == 1) {
		call printf ("replic = %d\n")
		    call pargi (NL_REPLI(nl))
	    } else {
		NL_REPLI(nl) = ival
	    }

	case RESTART: # :restart - List or set # of restarts.
	    call gargi (ival)
	    if (nscan() == 1) {
		call printf ("restart = %d\n")
		    call pargi (NL_RESTART(nl))
	    } else {
		NL_RESTART(nl) = ival
	    }

	case METHOD: # :method - List or set the minimization method.
	    call gargwrd (cmd, SZ_LINE)
	    if (nscan() == 1) {
		call printf ("method = %s\n")
		    call ic_gstr (ic, "algorithm", cmd, SZ_LINE)
		    call pargstr (cmd)
	    } else {
	        newm = strdic (cmd, cmd, SZ_LINE, METHODS)
	        if (newm != 0) {
	            call ic_pstr (ic, "algorithm", cmd)
	            call nl_puti (nl, "algorithm", newm)
	        } else 
	            call eprintf ("Non-recognizable method.\n")
	    }

	case VERBOSE: # :verbose - Set or list amoeba info control
	    call gargb (bval)
	    if (nscan() == 1) {
		call printf ("verbose = %b\n")
		    call pargb (nl_statb (nl, "verbose"))
	    } else
		call nl_putb (nl, "verbose", bval)

	case NGAUSS: # :ngauss - Set or list the number of Gaussians
	    ival = nl_stati (nl, "fitfunc")
	    if ((ival == GAUSSIANS) || (ival == CGAUSS)) {
	        call gargi (ival)
	        if (nscan() == 1) {
		    call printf ("ngauss = %d\n")
		        call pargi ((nl_stati (nl, "npar") - 2) / 3)
	        } else {
	            if (ival > 0) {
	                call icg_col4 (ic, nl, x, y, wts, 0,  ival * 3)
		        IC_NEWFUNCTION(ic) = YES
	            }
	        }
	    }

	case PCOMP: # :pcomp - Plot individual function components.
	    call gargb (bval)
	    if (nscan() == 1) {
		call printf ("pcomp = %b\n")
		if (IC_COMP(ic) == YES)
		    call pargb (YES)
	        else
		    call pargb (NO)
	    } else {
	        if (bval)
	             IC_COMP(ic) = YES
	        else
	             IC_COMP(ic) = NO
	        newgraph = YES
	    }

	default:
	    call eprintf ("Non-recognizable colon command.\n")
	}
end


# ICG_COL1 -- Change value of coefficient.

procedure icg_col1 (nl, index, val)

pointer	nl		# NLFIT pointer
int	index		# coefficient index
real	val		# value to assign to coefficient

int	nl_stati(), npar
pointer	coeff, sp

begin
	npar = nl_stati (nl, "npar")
	if ((index < 1) || (index > npar))
	    return

	# Get coefficients from NLFIT structure
	call smark (sp)
	call salloc (coeff, npar, TY_REAL)
	call nl_gcoeff (nl, Memr[coeff], npar)

	# Change desired coefficient
	Memr[coeff+index-1] = val

	# Save coefficients in NLFIT structure
	call nl_scoeff (nl, Memr[coeff], npar)
	call sfree (sp)
end


# ICG_COL2-- Get values associated with coefficient

int procedure icg_col2 (nl, index, val, err, flag)

pointer	nl		# NLFIT pointer
int	index		# coefficient index
real	val		# coefficient value
real	err		# coefficient error
bool	flag		# coefficient flag

int	nl_stati(), npar
pointer	coeff, errors, flags, sp

begin
	npar = nl_stati (nl, "npar")
	if ((index < 1) || (index > npar))
	    return (ERR)

	# Get coefficients from NLFIT structure
	call smark (sp)
	call salloc (coeff,  npar, TY_REAL)
	call salloc (errors, npar, TY_REAL)
	call salloc (flags, npar, TY_BOOL)
	call nl_gcoeff  (nl, Memr[coeff], npar)
	call nl_gerrors (nl, Memr[errors], npar)
	call nl_gflags  (nl, Memb[flags], npar)

	# Return values
	val  = Memr[coeff+index-1]
	err  = Memr[errors+index-1]
	flag = Memb[flags+index-1]

	call sfree (sp)
	return (OK)
end


# ICG_COL3 -- Change value of coefficient flag(s)

procedure icg_col3 (nl, index, bval)

pointer	nl		# NLFIT pointer
int	index		# coefficient index
bool	bval		# value to assign to coefficient flag

int	nl_stati(), npar, i
pointer	flags, sp

begin
	npar = nl_stati (nl, "npar")
	if ((index < 0) || (index > npar))
	    return

	# Get flags from NLFIT structure
	call smark (sp)
	call salloc (flags, npar, TY_BOOL)
	call nl_gflags (nl, Memb[flags], npar)

	# Change desired flag, or fix all flags.
	if (index == 0) {
	    do i = 0, npar-1
	        Memb[flags+i] = bval
	} else
	    Memb[flags+index-1] = bval

	# Save coefficients in NLFIT structure
	call nl_sflags (nl, Memb[flags], npar)
	call sfree (sp)
end


# ICG_COL4 -- Redefine fitting function in the nlfit package.

procedure icg_col4 (ic, nl, x, y, w, new1, nparg)

pointer	ic		# ICFIT pointer
pointer	nl		# NLFIT pointer
real	x[ARB]		# data arrays
real	y[ARB]
real	w[ARB]
int	new1		# new function
int	nparg		# new number of Gaussian coefficients

pointer	coeff, cflags
int	newfunc, oldfunc, oldnpar, npar, npts, i
real	amp

int	nl_stati()

begin
	newfunc = new1

	# Get old function data

	oldfunc = nl_stati (nl, "fitfunc")
	oldnpar = nl_stati (nl, "npar")
	npts    = nl_stati (nl, "npts")
	call malloc (coeff,  oldnpar, TY_REAL)
	call malloc (cflags, oldnpar, TY_BOOL)
	call nl_gcoeff (nl, Memr[coeff], oldnpar)
	call nl_gflags (nl, Memb[cflags], oldnpar)

	# Obtain estimate for amplitude of new component

	i = npts / 2
	while (w[i] == 0.)
	    i = i + 1
	amp = y[i]

	# Create new coefficient array. The action taken depends both on the 
	# old and new function types. The goal is to keep old coefficients, 
	# if possible. If not, put in their place reasonable values.
	# Some function type changes are forbidden, e.g. from/to Gaussians or
	# cgauss to/from other types.
	# To be more efficient, the following code is dependent on the way
	# the function coefficients are mapped on the array, for each 
	# specific functional form.

	switch (newfunc) {

	case POWERLAW:
	    switch (oldfunc) {
	    case BBODY, TWOBBODY:
	        npar = 3
	        call ic_pstr (ic, "function", "powerlaw")
	        call ic_puti (ic, "npar", npar)
	        Memr[coeff  + NL_PINDEX] = 1.
	        Memr[coeff  + NL_PAMPL]  = amp
	        # Keep same reference
	        Memb[cflags + NL_PINDEX] = true
	        Memb[cflags + NL_PAMPL]  = true
	        Memb[cflags + NL_PREF]   = false
	    case COMPOSITE:
	        npar = 3
	        call ic_pstr (ic, "function", "powerlaw")
	        call ic_puti (ic, "npar", npar)
	    default:
	        newfunc = oldfunc
	        npar    = oldnpar
	    } 

	case BBODY:
	    switch (oldfunc) {
	    case POWERLAW:
	        npar = 3
	        call ic_pstr (ic, "function", "bbody")
	        Memr[coeff  + NL_BTEMP] = 5000.
	        Memr[coeff  + NL_BAMPL]  = amp
		# Keep same reference
	        Memb[cflags + NL_BTEMP] = true
	        Memb[cflags + NL_BAMPL]  = true
	        Memb[cflags + NL_BREF]  = false
	    case COMPOSITE:
	        npar = 3
	        call ic_pstr (ic, "function", "bbody")
	        call ic_puti (ic, "npar", npar)
	        call amovr (Memr[coeff+npar], Memr[coeff], npar)
	        call amovi (Memb[cflags+npar], Memb[cflags], npar)
	    case TWOBBODY:
	        npar = 3
	        call ic_pstr (ic, "function", "bbody")
	        call ic_puti (ic, "npar", npar)
	    default:
	        newfunc = oldfunc
	        npar    = oldnpar
	    }

	case COMPOSITE:
	    switch (oldfunc) {
	    case POWERLAW:
	        npar = 6
	        call ic_pstr (ic, "function", "composite")
	        call ic_puti (ic, "npar", npar)
	        call realloc (coeff,  npar, TY_REAL)
	        call realloc (cflags, npar, TY_BOOL)
	        Memr[coeff  + NL_CTEMP]  = 5000.
	        Memr[coeff  + NL_CBAMPL] = amp
	        Memr[coeff  + NL_CBREF]  = Memr[coeff + NL_CPREF]
	        Memb[cflags + NL_CTEMP]  = true
	        Memb[cflags + NL_CBAMPL] = true
	        Memb[cflags + NL_CBREF]  = false
	    case BBODY:
	        npar = 6
	        call ic_pstr (ic, "function", "composite")
	        call ic_puti (ic, "npar", npar)
	        call realloc (coeff,  npar, TY_REAL)
	        call realloc (cflags, npar, TY_BOOL)
	        call amovr (Memr[coeff], Memr[coeff+3], 3)
	        call amovi (Memb[cflags], Memb[cflags+3], 3)
	        Memr[coeff  + NL_CINDEX] = 1.
	        Memr[coeff  + NL_CPAMPL] = amp
	        Memb[cflags + NL_CINDEX] = true
	        Memb[cflags + NL_CPAMPL] = true
	        Memb[cflags + NL_CPREF]  = false
	    case TWOBBODY:
	        npar = 6
	        call ic_pstr (ic, "function", "composite")
	        Memr[coeff  + NL_CINDEX] = 1.
	        Memr[coeff  + NL_CPAMPL] = amp
	        # Keep same reference
	        Memb[cflags + NL_CINDEX] = true
	        Memb[cflags + NL_CPAMPL] = true
	        Memb[cflags + NL_CPREF]  = false
	    default:
	        newfunc = oldfunc
	        npar    = oldnpar
	    }

	case TWOBBODY:
	    switch (oldfunc) {
	    case POWERLAW:
	        npar = 6
	        call ic_pstr (ic, "function", "twobbody")
	        call ic_puti (ic, "npar", npar)
	        call realloc (coeff,  npar, TY_REAL)
	        call realloc (cflags, npar, TY_BOOL)
	        Memr[coeff  + NL_TTEMP1] = 6000.
		Memr[coeff  + NL_TAMPL1] = amp
	        Memr[coeff  + NL_TTEMP2] = 5000.
	        Memr[coeff  + NL_TAMPL2] = amp
	        Memr[coeff  + NL_TREF2]  = Memr[coeff + NL_PREF]
	        Memb[cflags + NL_TTEMP1] = true
	        Memb[cflags + NL_TAMPL1] = true
	        Memb[cflags + NL_TREF1]  = false
	        Memb[cflags + NL_TTEMP2] = true
	        Memb[cflags + NL_TAMPL2] = true
	        Memb[cflags + NL_TREF2]  = false
	    case BBODY:
	        npar = 6
	        call ic_pstr (ic, "function", "twobbody")
	        call ic_puti (ic, "npar", npar)
	        call realloc (coeff,  npar, TY_REAL)
	        call realloc (cflags, npar, TY_BOOL)
	        call amovr (Memr[coeff], Memr[coeff+NL_TTEMP2], 3)
	        call amovi (Memb[cflags], Memb[cflags+NL_TTEMP2], 3)
	    case COMPOSITE:
	        npar = 6
	        call ic_pstr (ic, "function", "twobbody")
	        Memr[coeff  + NL_TTEMP1] = 5000.
	        Memr[coeff  + NL_TAMPL1] = amp
	        # Keep same reference
	        Memb[cflags + NL_TTEMP1] = true
	        Memb[cflags + NL_TAMPL1] = true
	        Memb[cflags + NL_TREF1]  = false
	    default:
	        newfunc = oldfunc
	        npar    = oldnpar
	    }

	case GAUSSIANS:
	    switch (oldfunc) {
	    case CGAUSS:
	        npar = oldnpar
	        call ic_pstr (ic, "function", "Gaussians")
	        if (oldnpar > 5) {
	            do i = 2, (oldnpar-2) / 3 {
	                Memr[coeff + NL_GAMPL(i)] = Memr[coeff + NL_GAMPL(i)] * Memr[coeff + NL_GAMPL(1)]
	                Memr[coeff + NL_GCENT(i)] = Memr[coeff + NL_GCENT(i)] + Memr[coeff + NL_GCENT(1)]
	            }
	        }
	    default:
	        newfunc = oldfunc
	        npar    = oldnpar
	    }

	case CGAUSS:
	    switch (oldfunc) {
	    case GAUSSIANS:
	        npar = oldnpar
	        call ic_pstr (ic, "function", "cgauss")
	        if (oldnpar > 5) {
	            do i = 2, (oldnpar-2) / 3 {
	                Memr[coeff + NL_GAMPL(i)] = Memr[coeff + NL_GAMPL(i)] / Memr[coeff + NL_GAMPL(1)]
	                Memr[coeff + NL_GCENT(i)] = Memr[coeff + NL_GCENT(i)] - Memr[coeff + NL_GCENT(1)]
	            }
	        }
	    default:
	        newfunc = oldfunc
	        npar    = oldnpar
	    }

	case 0:				# just changes number of Gaussians
	    switch (oldfunc) {
	    case GAUSSIANS, CGAUSS:
	        newfunc = oldfunc
	        npar = nparg + 2
	        call ic_puti (ic, "npar", npar)
	        if (oldnpar != npar) {
	            call realloc (coeff,  npar, TY_REAL)
	            call realloc (cflags, npar, TY_BOOL)
	        }
	        if (oldnpar < npar) {
	            call amovkr (0., Memr[coeff+oldnpar],  npar - oldnpar)
	            do i = oldnpar, npar-1
	                Memb[cflags+i] = false
	            call printf ("Don't forget to supply first guesses !")
	        }
	    default:
	        newfunc = oldfunc
	        npar    = oldnpar
	    }

	case USER,GALPROF:
	    newfunc = oldfunc
	    npar    = oldnpar

	default:
	    call error (0, "Non recognized function.")

	}

	call nl_init (nl, newfunc, Memr[coeff], Memb[cflags], npar, npts)

	call mfree (coeff,  TY_REAL)
	call mfree (cflags, TY_BOOL)
end

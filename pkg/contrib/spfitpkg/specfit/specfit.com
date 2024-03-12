# Common task parameters for SPECFIT

#Variables for the input data and errors
char	specname[SZ_LINE] # Name of current spectrum
int	npts		# Number of points in spectrum
real	itime		# Integration time
double	jd		# Julian Date
real	v0, v1		# Use to calculate variances for 1D images
pointer	spectrum	# Real pointer to location of data
pointer	lambdas		# Real pointer to location of wavelengths
pointer	errors		# Real pointer to location of associated errors

common /spec/ specname, npts, itime, jd, v0, v1, spectrum, lambdas, errors

#Variables defining the fitted sample
int	nfitpts		# Number of points used in the fit
int	nsamp		# Number of sample regions defined
real	sample[2, MAXSAMPLE] # Array defining start and end points of data to
			     # be used in the fit
pointer	infit		# Integer pointer to array of in/out of fit flags
pointer	fitsp		# Real pointer to fit values

common /samp/ nfitpts, nsamp, sample, infit, fitsp

#Variables describing the individual fit components

char	compkeys[SZ_CNAME,NKEYS] # Ascii keywords for component types
int	ncpar[NKEYS]		# Number of parameters for this component type
int	ncomp			# Number of components
int	comtyp[NCOMP]		# Type of component
int	parptr[PARPERC,NCOMP]	# Pointer to the parameters for each component
char 	comments[SZ_LINE,NCOMP] # Comments user puts in database file

common /comps/ ncpar, ncomp, comtyp, parptr, compkeys, comments

# Variables describing the fit parameters

int	npar		# Number of parameters in use
int	iptr[MAXPAR]	# Pointer to original locations of variables in fpar
int	ifix[MAXPAR]	# Flag for fixing (1) or freeing (0) a parameter
real	par0[MAXPAR]	# Parameter values
real	blim[MAXPAR]	# Hard lower limit for parameter
real	tlim[MAXPAR]	# Hard upper limit to parameter
real	step[MAXPAR]	# Step size for search
real	ptol[MAXPAR]	# Maximum allowed error in best fit
real	sigpar[MAXPAR,MAXPAR]	# One sigma error on the param.
				# Defined only for ifix=0

common /param/ npar, par0, blim, tlim, step, ptol, sigpar, ifix, iptr

# Variables for controlling the IMSL minimization routine

int	itr		# Maximum number of iterations
real	tolerance	# Fractional tolerance of the fit

common /contrl/ itr, tolerance

# Parameters that govern user interaction

bool	err_from_model
bool	interact
bool	debug
int	nlogfd
int	logfd[MAXLOG]

common /users/ err_from_model, interact, debug, nlogfd, logfd

# Variables describing the flux intervals of interest
char	linename[SZ_CNAME, NCOMP]	# Descriptive name for log file
int	ncont			# Component number for continuum
int	nlines			# Number of lines described
real	lam1[NCOMP]		# Starting wavelength of flux interval
real	lam2[NCOMP]		# Ending wavelength of flux interval
int	npieces[NCOMP]		# Number of components in each line
int	cmpnum[NCOMP, NCOMP]	# Fit component numbers for each

common /ansr/ linename, ncont, nlines, lam1, lam2, npieces, cmpnum

# Variables associated with grids of user line profiles
char	pmodname[SZ_FNAME,NMODMAX]
int	pfirst, pnmod, ppts
real	pmodkey[NMODMAX]
pointer	pw, pf

common /pmodel/ pfirst, pnmod, ppts, pmodname, pmodkey, pw, pf

# Variables associated with grids of user absorption models
char	amodname[SZ_FNAME,NMODMAX]
int	afirst, anmod, apts
real	amodkey[NMODMAX]
pointer	aw, af

common /amodel/ afirst, anmod, apts, amodname, amodkey, aw, af

# Variables associated with grids of user continuum models
char	cmodname[SZ_FNAME,NMODMAX]
int	cfirst, cnmod, cpts
real	cmodkey[NMODMAX]
pointer	cw, cf

common /cmodel/ cfirst, cnmod, cpts, cmodname, cmodkey, cw, cf

 
# Variables for use with the interrupts
int old_interrupt
common /interrupt_error/ old_interrupt
 

# Variable that tells which componenets are being plotted
int inplot[NCOMP]

common /list_of_plots_values/ inplot









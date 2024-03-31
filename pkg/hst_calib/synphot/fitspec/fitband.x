#* HISTORY *
#* B. Simon	28-Jul-94	original

# FITBAND -- Compute a least squares fit of a model to a bandpass

procedure fitband ()

#--
pointer	input		# observed bandpass
pointer	obsmode		# model bandpass
pointer	output		# table containing fitted bandpass
double	ftol		# fractional tolerance termination condition
int	maxiter		# maximum number of iterations
int	nprint		# number of iterations between diagnostic prints
bool	slow		# use slow method (simplex) to compute fit?
bool	equal		# use equal weighting on the data points?
double	var[9]		# free variables in least squares fit
pointer	grftable	# instrument graph table
pointer	cmptable	# component name table
real	hstarea		# total telescope area

int	nvar
pointer sp

bool	clgetb()
double	clgetd()
int	clgeti()
real	clgetr()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (input, SZ_FNAME, TY_CHAR)
	call salloc (obsmode, SZ_FNAME, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (grftable, SZ_FNAME, TY_CHAR)
	call salloc (cmptable, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgstr ("input", Memc[input], SZ_FNAME)
	call clgstr ("obsmode", Memc[obsmode], SZ_FNAME)
	call clgnone ("output", Memc[output], SZ_FNAME)

	ftol = clgetd ("ftol")
	maxiter = clgeti ("maxiter")
	nprint = clgeti ("nprint")
	slow = clgetb ("slow")
	equal = clgetb ("equal")
	call getfitvar (var, nvar)

	call clgstr ("grtbl", Memc[grftable], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptable], SZ_FNAME)
	hstarea = clgetr ("area") 
	call put_hstarea (hstarea)

	# Compute the values of the variables which minimize 
	# the least squares error of the fit

	call fitexpr (Memc[input], Memc[obsmode], Memc[output], ftol, 
		      maxiter, nprint, slow, equal, Memc[grftable], 
		      Memc[cmptable], var, nvar)

	# Write results to parameter file

	call putfitvar (var, nvar)
	call sfree (sp)
end

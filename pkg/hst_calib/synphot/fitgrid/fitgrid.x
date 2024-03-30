#* HISTORY *
#* B. Simon	21-Sep-94	original

# FITGRID -- Find the best fit to a spectrum from a grid of spectra

procedure fitgrid ()

#--
pointer	input		# observed spectrum 
pointer	spectrum	# list of spectra
pointer	output		# table containing fitted spectrum
pointer	vzero		# variable list
double	ftol		# fractional tolerance termination condition
int	maxiter		# maximum number of iterations
int	nprint		# number of iterations between diagnostic prints
bool	slow		# use slow method (simplex) to compute fit?
bool	equal		# use equal weighting on the data points?
pointer	grftable	# instrument graph table
pointer	cmptable	# component name table
real	hstarea		# total telescope area

pointer sp

bool	clgetb()
double	clgetd()
int	clgeti()
real	clgetr()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (input, SZ_FNAME, TY_CHAR)
	call salloc (spectrum, SZ_FNAME, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (vzero, SZ_FNAME, TY_CHAR)
	call salloc (grftable, SZ_FNAME, TY_CHAR)
	call salloc (cmptable, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgstr ("input", Memc[input], SZ_FNAME)
	call clgstr ("spectrum", Memc[spectrum], SZ_FNAME)
	call clgnone ("output", Memc[output], SZ_FNAME)
	call clgstr ("vzero", Memc[vzero], SZ_FNAME)

	ftol = clgetd ("ftol")
	maxiter = clgeti ("maxiter")
	nprint = clgeti ("nprint")
	slow = clgetb ("slow")
	equal = clgetb ("equal")

	call clgstr ("grtbl", Memc[grftable], SZ_FNAME)
	call clgstr ("cmptbl", Memc[cmptable], SZ_FNAME)
	hstarea = clgetr ("area") 
	call put_hstarea (hstarea)

	# Find the best fit to the input spectrum 
	# from all the spectra in the grid

	call bestgrid (Memc[input], Memc[spectrum], Memc[output], Memc[vzero],
		       ftol, maxiter, nprint, slow, equal, Memc[grftable], 
		       Memc[cmptable])


	call sfree (sp)
end

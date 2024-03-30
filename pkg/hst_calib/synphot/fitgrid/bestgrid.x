include	<tbset.h>
define	LENVAR		2

#* HISTORY *
#* B. Simon	21-Sep-94	original

# BESTGRID -- Find the best match to observed spectra from a grid of spectra

procedure bestgrid (input, spectrum, output, vzero, ftol, maxiter, 
		    nprint, slow, equal, grftable, cmptable)

char	input[ARB]	# i: file containing observed data
char	spectrum[ARB]	# i: synphot expression used as model
char	output[ARB]	# i: table containing intermediate results
char	vzero[ARB]	# i: list of values for variable zero
double	ftol		# i: chi squared termination condition
int	maxiter		# i: maximum number of iterations
int	nprint		# i: number of iterations between diagnostic prints
bool	slow		# i: use slow method (simplex) to find lsq fit?
bool	equal		# i: use equal weighting on the data points?
char	grftable[ARB]	# i: instrument graph table
char	cmptable[ARB]	# i: component name table
#--
double	var[LENVAR]
int	fd, ngrid, igrid, nwave, nvar, indegree, calcdegree, better, best
real	v0
pointer	sp, spec, spec2, units, pcode, chi2, scale
pointer	flux, sptr, wave, data, weight

string	nochoice  "Spectrum must contain at least two choices"
string	baddegree "Degree of input and spectrum do not match"

int	numlist(), numvzero(), nxtlist(), nxtvzero(), locpr()
pointer	rdlist(), open()
extern	getsynvar

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (spec, SZ_LINE, TY_CHAR)
	call salloc (spec2, SZ_LINE, TY_CHAR)
	call salloc (units, SZ_COLNAME, TY_CHAR)
	call salloc (pcode, SZ_LINE, TY_INT)

	# Initialize table cache and synphot variables

	call inisyntab
	call undefsynvar

	# Process the spectrum and vzero strings

	sptr = rdlist (spectrum)
	call rdvzero (vzero)

	# Get the total number of spectra

	ngrid = numlist (sptr) * numvzero ()
	if (ngrid < 2)
	    call printerr_str (nochoice, spectrum)

	call salloc (chi2, ngrid, TY_REAL)
	call salloc (scale, ngrid, TY_REAL)

	# Read the data from the input table

	call rdfittab (input, equal, indegree, Memc[units], 
		       nwave, wave, data, weight)

	call salloc (flux, nwave, TY_REAL)

	# Loop over each spectrum in the grid

	igrid = 0
	fd = open ("speclist", READ_WRITE, SPOOL_FILE)

	while (nxtlist (sptr, Memc[spec], SZ_LINE) != EOF) {

	    # Compile expression into pseudocode

	    call expcompile (Memc[spec], Memi[pcode], SZ_LINE)

	    # Loop over each value of vzero

	    while (nxtvzero (v0) != EOF) {

		# Write grid spectrum to spool file for later reference

		call fillexpr (Memc[spec], Memc[spec2], SZ_LINE)
		call fprintf (fd, "%s\n")
		call pargstr (Memc[spec2])

		# Compute flux of grid spectrum

		call syncalc (Memi[pcode], SZ_LINE, locpr(getsynvar), 
			      nwave, Memr[wave], grftable, cmptable, 
			      Memr[flux], calcdegree)

		if (indegree != calcdegree)
		    call printerr_str (baddegree, Memc[spec])

		# Compute chi squared and scaling factor for spectrum

		call specscale (nwave, Memr[wave], Memr[data], Memr[weight], 
				Memr[flux], Memr[scale+igrid], 
				Memr[chi2+igrid])
				
				

		igrid = igrid + 1
	    }
	}


	# Find the two spectra with the smallest chi squared values

	call rankgrid (ngrid, Memr[chi2], better, best)

	# Send diagnostic prints to user

	if (nprint > 0)
	    call reportgrid (fd, better, best, ngrid, Memr[chi2], Memr[scale])

	# Do a least squares fit between these two spectra

	call selectgrid (fd, better, best, Memr[scale], 
			nvar, var, Memc[spec], SZ_LINE)

	call tweengrid (Memc[spec], nwave, Memr[wave], Memr[data],
			Memr[weight], output, indegree, Memc[units], 
			ftol, maxiter, nprint, slow, grftable, cmptable, 
			var, nvar)

	# Close files and release memory

	call close (fd)
	call clssyntab

	call mfree (wave, TY_REAL)
	call mfree (data, TY_REAL)
	call mfree (weight, TY_REAL)

	call freelist (sptr)
	call sfree (sp)
end

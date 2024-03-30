include	<tbset.h>

#* HISTORY *
#* B. Simon	29-Jul-94	original

# FITFUNC -- Coumpute a function of the fit variables

procedure fitfunc (ndata, nvar, var, resid, iflag)

char	input[1]	# i: file containing observed data
char	expr[1]		# i: synphot expression used as model
char	output[1]	# i: output table name
char	grftable[1]	# i: instrument graph table
char	cmptable[1]	# i: component name table
bool	equal		# i: use equal weighting on the observed data?
int	ndata		# i: number of data points
int	nvar		# i: number of fit variables
double	var[ARB]	# i: fit variables
double	resid[ARB]	# o: residual of the fit
int	iflag		# i: print and error flag
#--
int	irep		# current repetition of the call to fitfunc
double	chisq		# sum of squares of the residual to the fit
int	indegree	# degree of the data in the input file
int	nwave		# length of wavelength array
pointer	fluxunits	# flux units of input table
pointer	grftab		# graph table name
pointer	cmptab		# component lookup table name
pointer	wave		# wavelengths at which fit is evaluated
pointer	data		# observed data to be fit
pointer	weight		# weights used in computing residuals
pointer	command		# synphot command containing the expression to fit
pointer	pcode		# pseudocode use to evaluate the expression
pointer	model		# result of the expression after evaluation

int	ic, jc, ivar, calcdegree, ncol
pointer	fd, sp, command2, line
real	value

string	diagfmt   "irep = %d chisq = %8g exp = %s\n"
string	finalfmt  "\nFinal solution:\n"
string	badunits  "Units of input data and expression do not agree"

extern	getsynvar, putfile
int	getline(), envgeti(), gstrcpy()
pointer	open(), locpr()

begin
	# Check value of print flag, print intermediate results
	# or calculate fitted function and its residuals

	if (iflag == 0) {
	    call smark (sp)
	    call salloc (command2, SZ_FNAME, TY_CHAR)
	    call fillexpr (Memc[command], Memc[command2], SZ_FNAME)

	    # Write intermediate results to STDERR

	    call fprintf (STDERR, diagfmt)
	    call pargi (irep)
	    call pargd (chisq)
	    call pargstr (Memc[command2])

	    call sfree (sp)

        } else if (iflag > 0) {
	    # Set values of variables in fit

	    if (iflag == 1)
		irep = irep + 1

	    do ivar = 1, nvar {
		value = var[ivar]
		call putsynvar (ivar, value)
	    }

	    # Compute fitted function

	    call syncalc (Memi[pcode], SZ_COMMAND, locpr(getsynvar), 
			  ndata, Memr[wave], Memc[grftab], Memc[cmptab], 
			  Memr[model], calcdegree)

	    if (indegree != calcdegree)
		call printerr_str (badunits, "fitfunc")

	    # Calculate residuals from fit

	    call getresid (ndata, Memr[model], Memr[data], Memr[weight], 
			   resid, chisq)
        }

        return

	# SETFITFUNC -- Initialize variables used in the fit

        entry setfitfunc (input, expr, grftable, cmptable, 
			  equal, ndata)

	call smark (sp)
	call salloc (line, SZ_LINE, TY_CHAR)

	# Initialize global variables used in the fit

	call malloc (fluxunits, SZ_COLUNITS, TY_CHAR)
	call malloc (grftab, SZ_FNAME, TY_CHAR)
	call malloc (cmptab, SZ_FNAME, TY_CHAR)
	call malloc (command, SZ_COMMAND, TY_CHAR)
	call malloc (pcode, SZ_COMMAND, TY_INT)

	# Read the observed data from the input file

	call rdfittab (input, equal, indegree, Memc[fluxunits], 
		       ndata, wave, data, weight)

	call malloc (model, ndata, TY_REAL)
	nwave = ndata

	irep = 0
	chisq = 0.0
	call lastfile (grftable, Memc[grftab], SZ_FNAME)
	call lastfile (cmptable, Memc[cmptab], SZ_FNAME)

	# Search for first non-white character in expression

	for (ic = 1; expr[ic] != EOS; ic = ic + 1)
	    if (expr[ic] > ' ')
		break

	# Copy expression into command, possibly reading from file

	jc = 0
	if (expr[ic] != '@') {
	    call strcpy (expr[ic], Memc[command], SZ_COMMAND)

	} else {
	    fd = open (expr[ic+1], READ_ONLY, TEXT_FILE)

	    while (getline (fd, Memc[line]) != EOF) {
		jc = jc + gstrcpy (Memc[line], Memc[command+jc], SZ_COMMAND-jc)
		Memc[command+jc-1] = ' '
		Memc[command+jc] = ' '
		jc = jc + 1
	    }

	    Memc[command+jc] = EOS
	    call close (fd)
	}

	# Compile the command into pseudocode

	call expcompile (Memc[command], Memi[pcode], SZ_COMMAND)

	call sfree (sp)
	return


	# CLSFITFUNC -- Clean up after fitfunc

	entry clsfitfunc (output, nvar, var)

	call smark (sp)
	call salloc (command2, SZ_COMMAND, TY_CHAR)

	# Fill the expression with the final values of the fit variables
	# and write the results to STDERR

	do ivar = 1, nvar {
	    value = var[ivar]
	    call putsynvar (ivar, value)
	}

	call fillexpr (Memc[command], Memc[command2], SZ_COMMAND)

	iferr {
	    ncol = envgeti ("ttyncols") - 1
	} then {
	    ncol = 79
	}

	call fprintf (STDERR, finalfmt)
	call wrtcommand (putfile, STDERR, ncol, Memc[command2])

	# Write the fitted function to the output table

	call wrtfittab (output, Memc[grftab], Memc[cmptab], Memc[fluxunits], 
			Memc[command2], nvar, var, nwave, Memr[wave], 
			SZ_COMMAND, Memi[pcode])

	# Free memory used to hold global variables

	call mfree (fluxunits, TY_CHAR)
	call mfree (grftab, TY_CHAR)
	call mfree (cmptab, TY_CHAR)
	call mfree (wave, TY_REAL)
	call mfree (data, TY_REAL)
	call mfree (weight, TY_REAL)
	call mfree (command, TY_CHAR)
	call mfree (pcode, TY_INT)
	call mfree (model, TY_REAL)

	call sfree (sp)

	return
end

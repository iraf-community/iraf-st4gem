#* HISTORY *
#* B.Simon	30-Sep-94	Original

# GRIDFUNC -- Evaluate function and compute residuals of grid spectra

procedure gridfunc (nwave, nvar, var, resid, iflag)

int	nwave		# i: number of data points
int	nvar		# i: number of fit variables
double	var[ARB]	# i: fit variables
double	resid[ARB]	# o: residual of the fit
int	iflag		# i: print and error flag
char 	spec[1]		# i: spectrum expression
char	grftable[1]	# i: instrument graph table
char	cmptable[1]	# i: component name table
real	wave[1]		# i: wavelengths at which fit is evaluated
real	data[1]		# i: observed data to be fit
real	weight[1]	# i: weights used in computing residuals
char	output[1]	# i: output file name
char	units[1]	# i: input (and output) units
#--
int	irep		# current repetition of the call to fitfunc
double	chisq		# sum of squares of the residual to the fit
pointer	grftab		# graph table name
pointer	cmptab		# component lookup table name
pointer	svwave		# wavelengths at which fit is evaluated
pointer	svdata		# observed data to be fit
pointer	svweight        # weights used in computing residuals
pointer	command		# synphot command containing the expression to fit
pointer	pcode		# pseudocode use to evaluate the expression

int	ivar, degree, ncol
pointer	sp, command2, model
real	value

string	diagfmt   "irep = %d chisq = %8g exp = %s\n"
string	finalfmt  "\nFinal solution:\n"

extern	getsynvar, putfile
int	envgeti()
pointer	locpr()

begin
	# Allocate temporary variables

	call smark (sp)
	call salloc (command2, SZ_FNAME, TY_CHAR)
	call salloc (model, nwave, TY_REAL)

	# Set values of variables in fit

	do ivar = 1, nvar {
	    value = var[ivar]
	    call putsynvar (ivar, value)
	}

	if (iflag == 0) {
	    # Write intermediate results to STDERR

	    call fillexpr (Memc[command], Memc[command2], SZ_FNAME)

	    call fprintf (STDERR, diagfmt)
	    call pargi (irep)
	    call pargd (chisq)
	    call pargstr (Memc[command2])

        } else if (iflag > 0) {
	    # Compute fitted function

	    if (iflag == 1)
		irep = irep + 1

	    call syncalc (Memi[pcode], SZ_COMMAND, locpr(getsynvar), 
			  nwave, Memr[svwave], Memc[grftab], Memc[cmptab], 
			  Memr[model], degree)

	    # Calculate residuals from fit

	    call getresid (nwave, Memr[model], Memr[svdata], Memr[svweight], 
			   resid, chisq)
	}

	call sfree (sp)
	return

	# SETGRIDFUNC -- Initialize variables used in the fit

        entry setgridfunc (spec, grftable, cmptable, nwave, wave, data, 
			  weight)

	# Initialize global variables used in the fit

	call malloc (grftab, SZ_FNAME, TY_CHAR)
	call malloc (cmptab, SZ_FNAME, TY_CHAR)
	call malloc (command, SZ_COMMAND, TY_CHAR)
	call malloc (pcode, SZ_COMMAND, TY_INT)

	call malloc (svwave, nwave, TY_REAL)
	call malloc (svdata, nwave, TY_REAL)
	call malloc (svweight, nwave, TY_REAL)

	irep = 0
	chisq = 0.0

	call amovr (wave, Memr[svwave], nwave)
	call amovr (data, Memr[svdata], nwave)
	call amovr (weight, Memr[svweight], nwave)

	call strcpy (spec, Memc[command], SZ_COMMAND)

	call lastfile (grftable, Memc[grftab], SZ_FNAME)
	call lastfile (cmptable, Memc[cmptab], SZ_FNAME)

	# Compile the command into pseudocode

	call expcompile (spec, Memi[pcode], SZ_COMMAND)
	return


	# CLSFITFUNC -- Clean up after fitgrid

	entry clsgridfunc (output, units, nvar, var, nwave, wave)

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

	call wrtfittab (output, Memc[grftab], Memc[cmptab], units, 
			Memc[command2], nvar, var, nwave, wave, 
			SZ_COMMAND, Memi[pcode])

	# Free memory used to hold global variables

	call mfree (grftab, TY_CHAR)
	call mfree (cmptab, TY_CHAR)
	call mfree (command, TY_CHAR)
	call mfree (pcode, TY_INT)
	call mfree (svwave, TY_REAL)
	call mfree (svdata, TY_REAL)
	call mfree (svweight, TY_REAL)

	call sfree (sp)
	return

end

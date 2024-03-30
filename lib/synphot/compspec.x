# COMPSPEC -- Spectrum expression evaluator

procedure compspec (command, iw, graphtab, comptab, nwave, wave, spec, form)

char	command[ARB]	# i: string containing expression to evaluate
int	iw		# i: position in string where expression starts
char	graphtab[ARB]	# i: graph table name
char	comptab[ARB]	# i: component lookup table name
int	nwave		# i: number of wavelengths
pointer	wave		# u: wavelength array
pointer	spec		# o: spectrum array, calling array must free
char	form[ARB]	# o: form (units) of spectrum
#--
int	units, intersect
pointer	sp, pcode, expr

data	intersect  / YES /
string	notspec   "Expression is not a spectrum"

errchk	syncompile, inisyntab, synwave, syncalc, synphoterr

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (pcode, SZ_COMMAND, TY_INT)

	# Allocate memory for output array

	call malloc (spec, nwave, TY_REAL)

	# Convert expression to pseudocode for calculator

	call syncompile (command[iw], Memi[pcode], SZ_COMMAND)

	# If the wavelength set is not supplied, compute it
	# (flagged by putting an INDEF in the first element)

	call inisyntab
	if (IS_INDEFR (Memr[wave])) {
	    call synwave (intersect, Memi[pcode], 1, SZ_COMMAND, 
			  graphtab, comptab, Memr[wave], nwave)
	}

	# Calculate spectrum from pseudocode

	call syncalc (Memi[pcode], SZ_COMMAND, NULL, nwave, Memr[wave], 
		      graphtab, comptab, Memr[spec], units)
	call clssyntab

	# Check output units. Not one means expression does not
	# evaluate into a spectrum

	if (units != 1) {
	    call salloc (expr, 40, TY_CHAR)
	    call strcpy (command, Memc[expr], 40)

	    call synphoterr (notspec, Memc[expr])
	}

	call strcpy ("photlam", form, 8)
	call sfree (sp)
end

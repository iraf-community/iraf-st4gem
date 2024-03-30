# COMPBAND -- Passband expression evaluator

procedure compband (command, iw, graphtab, comptab, nwave, wave, band)

char	command[ARB]	# i: string containing expression to evaluate
int	iw		# i: position in string where expression starts
char	graphtab[ARB]	# i: graph table name
char	comptab[ARB]	# i: component lookup table name
int	nwave		# i: number of wavelengths
pointer	wave		# u: wavelength array
pointer	band		# o: passband array, calling subroutine must free
#--
bool	logspace
int	units, intersect
pointer	sp, pcode, expr, error

data	intersect  / YES /
data	logspace   / true /

string	notband   "Not a passband"

bool	is_simple()
errchk	syncompile, inisyntab, synwave, syncalc, synphoterr

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (pcode, SZ_COMMAND, TY_INT)
	call salloc (expr, SZ_FNAME, TY_CHAR)

	# Allocate memory for output array

	call malloc (band, nwave, TY_REAL)

	if (is_simple (command[iw])) {
	    # Bypass calculator for simple bandpass expressions

	    call salloc (error, nwave, TY_REAL)

	    if (IS_INDEFR (Memr[wave])) {
		call getbandx (command[iw], graphtab, comptab, logspace, 
			       nwave, Memr[wave], Memr[band], Memr[error])
	    } else {
		call evalbandx (command[iw], nwave, Memr[wave], graphtab, 
				comptab, Memr[band], Memr[error])
	    }

	} else {
	    # Convert expression to pseudocode for calculator

	    call syncompile (command[iw], Memi[pcode], SZ_COMMAND)

	    # If the wavelength set is not supplied, compute it
	    # (flagged by putting an INDEF in the first element)

	    call inisyntab

	    if (IS_INDEFR (Memr[wave]))
		call synwave (intersect, Memi[pcode], 1, SZ_COMMAND, 
			      graphtab, comptab, Memr[wave], nwave)


	    # Calculate bandpass from pseudocode

	    call syncalc (Memi[pcode], SZ_COMMAND, NULL, nwave, Memr[wave],
			  graphtab, comptab, Memr[band], units)
	    call clssyntab

	    # Check output units. Non-zero means expression does not
	    # evaluate into a passband

	    if (units != 0) {
		call strcpy (command, Memc[expr], SZ_FNAME)
		call synphoterr (notband, Memc[expr])
	    }
	}

	call sfree (sp)
end

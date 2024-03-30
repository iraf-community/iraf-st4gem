include	<tbset.h>

#* HISTORY *
#* B.Simon	16-Sep-94	original

# WRTFITTAB -- Create the least squares fit output table

procedure wrtfittab (output, grftab, cmptab, fluxunits, command, 
		     nvar, var, nwave, wave, maxcode, pcode)

char	output[ARB]	# i: output table name
char	grftab[ARB]	# i: graph table name
char	cmptab[ARB]	# i: component lookup table name
char	fluxunits[ARB]	# i: flux units of output table
char	command[ARB]	# i: synphot command
int	nvar		# i: number of fit variables
double	var[ARB]	# i: fit variables
int	nwave		# i: length of wavelength and data arrays
real	wave[ARB]	# i: wavelengths at which fit is evaluated
int	maxcode		# i: length of pseudocode array
int	pcode[ARB]	# i: pseudocode used to evaluate the expression
#--
int	ivar, degree, done
pointer	sp, data, out, wv, dt
real	value

string	blank  " "

extern	getsynvar
extern	puthead
int	phottoany()
pointer	tbtopn(), locpr()

begin
	if (output[1] == EOS)
	    return

	call smark (sp)
	call salloc (data, nwave, TY_REAL)

	# Compute fitted function

	do ivar = 1, nvar {
	    value = var[ivar]
	    call putsynvar (ivar, value)
	}

	call syncalc (pcode, maxcode, locpr(getsynvar), nwave, wave, 
		      grftab, cmptab, Memr[data], degree)

	if (degree == 1)
	    done = phottoany (fluxunits, nwave, wave, Memr[data])

	# Open output table

	out = tbtopn (output, NEW_FILE, NULL)
	call tbpset (out, TBL_MAXPAR, 25)

        call tbcdef (out, wv, "WAVELENGTH", "angstroms", blank, TY_REAL, 1, 1)
	if (degree == 0) {
	    call tbcdef (out, dt, "THROUGHPUT", blank, blank, TY_REAL, 1, 1)
	} else {
	    call tbcdef (out, dt, "FLUX", fluxunits, blank, TY_REAL, 1, 1)
	}

	call tbtcre (out)

	call tbhadt (out, "grftable", grftab)
	call tbhadt (out, "cmptable", cmptab)

	call wrtcommand (puthead, out, SZ_FNAME, command)

	call tbcptr (out, wv, wave, 1, nwave)
	call tbcptr (out, dt, Memr[data], 1, nwave)

	call tbtclo (out)
	call sfree (sp)
	return

end

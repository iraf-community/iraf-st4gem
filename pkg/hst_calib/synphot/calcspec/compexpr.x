include <tbset.h>
include <synphot.h>
define	NPHOT		7

# COMPEXPR -- Compute the result of a synphot expression

procedure compexpr (command, output, outform, vzero, 
		    wavtable, grftable, cmptable)

char	command[ARB]	# i: Expression to calculate
char	output[ARB]	# i: Output table name
char	outform[ARB]	# i: Form of output spectrum
char	vzero[ARB]	# i: Variable list
char	wavtable[ARB]	# i: Table containing wavelength array
char	grftable[ARB]	# i: Instrument graph table
char	cmptable[ARB]	# i: Component name table
#--
int	nexp, iexp, nwave, degree, done, iphot
real	area, v0, lambda, tlambda, phot[NPHOT]
pointer	tp, sp, expr, expr2, keywrd, tbname, pcode
pointer	cmd, wave, result, wptr, cptr

string	nilstr  ""
string	isempty "Spectrum file is empty"
string	keylist "uresp,pivwv,bandw,tpeak,equvw,rectw,emflx"

extern	getsynvar
int	numlist(), numvzero(), nxtlist(), nxtvzero(), phottoany(), word_find()
pointer	rdlist(), tbtopn(), locpr()
real	funit(), pivlam(), rmslam(), peaklam()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (expr, SZ_LINE, TY_CHAR)
	call salloc (expr2, SZ_FNAME, TY_CHAR)
	call salloc (keywrd, SZ_FNAME, TY_CHAR)
	call salloc (tbname, SZ_FNAME, TY_CHAR)
	call salloc (pcode, SZ_LINE, TY_INT)

	# Initialize table cache and synphot variables

	call inisyntab
	call undefsynvar

	# Process the command and vzero strings

	cmd = rdlist (command)
	call rdvzero (vzero)

	# Get the total number of expressions

	nexp = numlist (cmd) * numvzero ()
	if (nexp <= 0)
	    call printerr_str (isempty, command)

	# Create the output table

	tp = tbtopn (output, NEW_FILE, NULL)
	call tbpset (tp, TBL_ROWLEN, nexp+1)
	call tbpset (tp, TBL_MAXCOLS, nexp+1)
	call tbpset (tp, TBL_MAXPAR, (NPHOT+1)*nexp+6)
	call tbcdef (tp, wptr, "WAVELENGTH", "angstroms", nilstr, 
		     TY_REAL, 1, 1)
	call tbtcre (tp)

	# Add header parameters

	call lastfile (grftable, Memc[tbname], SZ_FNAME)
	call tbhadt (tp, "grftable", Memc[tbname])

	call lastfile (cmptable, Memc[tbname], SZ_FNAME)
	call tbhadt (tp, "cmptable", Memc[tbname])

	iexp = 1
	# Loop over each expression

	while (nxtlist (cmd, Memc[expr], SZ_LINE) != EOF) {

	    # Compile expression into pseudocode

	    call expcompile (Memc[expr], Memi[pcode], SZ_LINE)

	    # Get wavelength set on first pass through loop

	    if (iexp == 1) {
		call getwavelen (wavtable, grftable, cmptable, Memi[pcode], 
				 1, SZ_LINE, wave, nwave)
		call tbcptr (tp, wptr, Memr[wave], 1, nwave)
		call salloc (result, nwave, TY_REAL)
	    }

	    # Loop over each value of vzero

	    while (nxtvzero (v0) != EOF) {

		# Calculate result of expression

		call syncalc (Memi[pcode], SZ_COMMAND, locpr(getsynvar), 
			      nwave, Memr[wave], grftable, cmptable, 
			      Memr[result], degree)

		# Set appropriate column name  and convert units

		if (degree == 0) {
		    call makename ("THROUGHPUT", iexp, nexp, 
				   Memc[tbname], SZ_FNAME)
		} else {
		    done = phottoany (outform, nwave, Memr[wave], Memr[result])
		    call makename ("FLUX", iexp, nexp, Memc[tbname], SZ_FNAME)
		}

		# Calculate the photometric paraneters from the throughput

		# Write the result

		call tbcdef (tp, cptr, Memc[tbname], outform, nilstr, 
			     TY_REAL, 1, 1)
		call tbcptr (tp, cptr, Memr[result], 1, nwave)

		# Write the header keywords

		call fillexpr (Memc[expr], Memc[expr2], SZ_FNAME)
		call makename ("expr", iexp, nexp, Memc[tbname], SZ_FNAME)
		call tbhadt (tp, Memc[tbname], Memc[expr2])

		if (degree == 0) {
		    call get_hstarea (area)
		    call tbhadr (tp, "aperarea", area)
		    call tbhadr (tp, "zeropt", STZERO)

		    phot[1] = funit (area, nwave, Memr[wave], Memr[result])
		    phot[2] = pivlam (nwave, Memr[wave], Memr[result])
		    phot[3] = rmslam (nwave, Memr[wave], Memr[result])
		    phot[4] = peaklam (nwave, Memr[result])

		    call widthlam (nwave, Memr[wave], Memr[result], 
				   phot[5], phot[6])

		    lambda = INDEFR
		    call monolam (nwave, Memr[wave], Memr[result], lambda,
				  tlambda, phot[7])

		    do iphot = 1, NPHOT {
			if (word_find (iphot, keylist, Memc[keywrd], 
				       SZ_FNAME) < 0)
			    break

			call makename (Memc[keywrd], iexp, nexp, 
				       Memc[tbname], SZ_FNAME)
			call tbhadr (tp, Memc[tbname], phot[iphot])
		    }
		}

		iexp = iexp + 1
	    }
	}

	# Close files and release memory

	call clssyntab
	call tbtclo (tp)

	call mfree (wave, TY_REAL)
	call freelist (cmd)
	call sfree (sp)
end

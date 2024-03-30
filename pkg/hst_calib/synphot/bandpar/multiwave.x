include	<tbset.h>
include	<synphot.h>

#* HISTORY *
#* B.Simon	07-Oct-94	original
#* B.Simon	21-Aug-00	revides to call synwave

# MULTIWAVE -- Compute a wavelength set from multiple expressions

procedure multiwave (wavtable, grftable, cmptable, list, wave, nwave)
			
char	wavtable[ARB]	# i: Table containing wavelength array
char	grftable[ARB]	# i: Instrument graph table
char	cmptable[ARB]	# i: Component name table
pointer	list		# i: list of expressions
pointer	wave		# o: wavelength array (caller must deallocate)
int	nwave		# o: number of wavelengths
#--
int	ncode, intersect
pointer	tp, sp, expr, pcode, ocode

data	intersect  / NO /

int	tbpsta(), numlist(), nxtlist()
pointer	tbtopn()

begin
	# Allocate memory for temporary strings

	ncode = numlist (list)

	call smark (sp)
	call salloc (expr, SZ_LINE, TY_CHAR)
	call salloc (pcode, ncode*SZ_LINE, TY_INT)

	if (wavtable[1] != EOS) {
	    # Wavelength table specified, read wavelength array from it

	    tp = tbtopn (wavtable, READ_ONLY, NULL)
	    nwave = tbpsta (tp, TBL_NROWS)
	    call malloc (wave, nwave, TY_REAL)

	    call rdwave (tp, nwave, Memr[wave])
	    call tbtclo (tp)

	} else {
	    # No wavelength table, calculate wavelength set from expressions

	    # Get each expression and calculate its wavelength range
	    # the total range is the union of the individual ranges

	    ocode = pcode
	    while (nxtlist (list, Memc[expr], SZ_LINE) != EOF) {
		call expcompile (Memc[expr], Memi[ocode], SZ_LINE)
		ocode = ocode + SZ_LINE
	    }

	    call synwave (intersect, Memi[pcode], ncode, SZ_LINE, 
			  grftable, cmptable, wave, nwave)
	}

	call sfree (sp)
end

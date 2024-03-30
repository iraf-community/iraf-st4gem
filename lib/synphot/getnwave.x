include	<tbset.h>
include "syncalc.h"
include "libsynphot.h"

# GETNWAVE -- Get maximum number of wavelengths from tables used in expression

int procedure getnwave (pcode, maxcode)

int	pcode[ARB]	# i: pseudocode used by calculator
int	maxcode		# i: maximum length of pseudocode
#--
int	istr, nwave, icode, nrow
pointer	tp, wv, sp, strbuf,tabname, thruname, str

string	name	 "getnwave"
string	wavecol  "WAVELENGTH"
string	badcode  "Illegal code in expression evaluator"

int	tbcigi(), tbpsta()
pointer	opnsyntab()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (strbuf, maxcode, TY_CHAR)
	call salloc (tabname, SZ_FNAME, TY_CHAR)
	call salloc (thruname, SZ_FNAME, TY_CHAR)

	istr = 0
	nwave = 0
	icode = 1

	# Skip over all tokens that are not filenames
	# Set number of wavelengths to max rows in any table

	while (pcode[icode] != Y_DONE) {
	    switch (pcode[icode]) {
	    case Y_STR, Y_NUM:
		icode = icode + 1
		call calcstr (pcode, icode, strbuf, istr, str)
		icode = icode + 1

	    case Y_FILE:
		icode = icode + 1

		# Get table name

		call calcstr (pcode, icode, strbuf, istr, str)
		icode = icode + 1

		call breakcomp (Memc[str], Memc[tabname], 
				Memc[thruname], SZ_FNAME)

		# Check for wavelength column. If it exists,
		# we have a spectrum or throughput table

		tp = opnsyntab (Memc[tabname])
		call syncolptr (tp, wavecol, 1, wv)

		if (tbcigi (wv, TBL_COL_DATATYPE) > 0) {
		    nrow = tbpsta (tp, TBL_NROWS)
		    nwave = max (nwave, nrow)
		}

	    case Y_VAR:
		icode = icode + 2

	    case Y_FUNC:
		icode = icode + 3

	    case Y_ADD, Y_SUB, Y_MUL, Y_DIV, Y_NEG:
		icode = icode + 1

	    default:
		call synphoterr (badcode, name)
	    }
	}

	nwave = max (LENWAVE, nwave)
	return (nwave)

end

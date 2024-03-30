include "syncalc.h"
include "libsynphot.h"

define	MAXSTACK	64

# CALCFILES -- Calculate filenames used in expression

procedure calcfiles (pcode, maxcode, graphtab, comptab, logfile)

int	pcode[ARB]	# i: pseudocode used by calculator
int	maxcode		# i: maximum length of pseudocode
char    graphtab[ARB]   # i: graph table name
char    comptab[ARB]    # i: component lookup table name
char	logfile[ARB]	# i: log file name
#--
int	ic, nc, top, bot, istr, icode
int	iarg, narg, ftype, log
pointer	sp, nilstr, temp, strbuf, str
pointer	valstk[MAXSTACK]
real	value

string	name	 "calcfiles"
string	catfile  CAT_TEMPLATE
string	badexpr  "Synthetic photometry expression too complex"
string	badcode  "Illegal code in expression evaluator"

int	gstrcpy(), open(), ctor()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (nilstr, 1, TY_CHAR)
	call salloc (temp, SZ_LINE, TY_CHAR)
	call salloc (strbuf, maxcode, TY_CHAR)

	Memc[nilstr] = EOS

	# Top is current top of stack

	top = 0
	istr = 0
	icode = 1
	log = open (logfile, APPEND, TEXT_FILE)

	while (pcode[icode] != Y_DONE) {

	    # Make room for results of current operation

	    top = top + 1
	    if (top > MAXSTACK)
		call synphoterr (badexpr, name)

	    switch (pcode[icode]) {
	    case Y_STR:
		# String stored in code.

		narg = 0
		icode = icode + 1
		call calcstr (pcode, icode, strbuf, istr, str)

		valstk[top] = str

	    case Y_FILE:
		# Filename stored in code.

		narg = 0
		icode = icode + 1
		call calcstr (pcode, icode, strbuf, istr, str)

		call fprintf (log, "#Filename in expression:\n%s\n\n")
		call pargstr (Memc[str])

		valstk[top] = str

	    case Y_NUM:
		# Number is stored in code as string. 

		narg = 0
		icode = icode + 1
		call calcstr (pcode, icode, strbuf, istr, str)

		valstk[top] = str

	    case Y_VAR:
		# Variable number follows in code

		narg = 0
		icode = icode + 1

		valstk[top] = nilstr

	    case Y_FUNC:
		# Function type and number of arguments stored in code.

		icode = icode + 1
		ftype = pcode[icode]
		icode = icode + 1
		narg = pcode[icode]

		valstk[top] = nilstr

		ic = 0
		do iarg = top-narg, top-1 {
		    ic = ic + gstrcpy (Memc[valstk[iarg]], Memc[temp+ic], 
				       SZ_LINE - ic)
		    Memc[temp+ic] = ','
		    ic = ic + 1
		}
		Memc[temp+ic] = EOS

		switch (ftype) {
		case FN_BAND:
		    call fprintf (log, "#Throughput table names:\n")
		    call bandname (Memc[temp], graphtab, comptab, log)
		    call fprintf (log, "\n")

		case FN_CAT:
		    call fprintf (log, "#Catalog table name:\n")
		    call catname (Memc[temp], catfile, NO, log)
		    call fprintf (log, "\n")

		case FN_GRID:
		    ic = 1
		    nc = ctor (Memc[valstk[top-1]], ic, value)

		    call fprintf (log, "#Grid table names:\n")
		    call gridname (Memc[valstk[top-2]], value, log)
		    call fprintf (log, "\n")

		case FN_ICAT:
		    call fprintf (log, "#Catalog table names:\n")
		    call catname (Memc[temp], catfile, YES, log)
		    call fprintf (log, "\n")

		case FN_SPEC:
		    call fprintf (log, "#Spectrum name:\n%s\n\n")
		    call pargstr (Memc[valstk[top-1]])

		case FN_THRU:
		    call fprintf (log, "#Throughput table name:\n%s\n\n")
		    call pargstr (Memc[valstk[top-1]])
		}

	    case Y_ADD, Y_SUB, Y_MUL, Y_DIV:
		# Four arithmetic functions

		narg = 2
		valstk[top] = nilstr

	    case Y_NEG:
		narg = 1
		valstk[top] = nilstr

	    default:
		call synphoterr (badcode, name)
	    }

	    if (narg > 0) {
 		bot = top - narg
		valstk[bot] = valstk[top]
		top = bot
	    }

	    icode = icode + 1
	}

	call close (log)
end



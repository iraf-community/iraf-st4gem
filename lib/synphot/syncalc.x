include	<tbset.h>
include	<error.h>
include "syncalc.h"
include	"libsynphot.h"

define	TY_STRING	1
define	TY_NUMBER	2
define	TY_VECTOR	3
define	MAXSTACK	64

# SYNCALC -- Calculate the result of a synthetic photometry expression

# This procedure takes an expression that has been compiled into RPN
# code and calculates the results of the expression. Filenames and other
# constants the expression operates on are embedded in the code. 
# Intermediate results of the calculation are stored on a stack. The
# stack tracks three pieces of information: a pointer to the memory 
# holding the intermediate result, the type of the intermediate result, 
# and its units. The units are expressed as an integer power of FLAM, 
# the form used internally for calculations.

procedure syncalc (pcode, maxcode, getvar, nwave, wave,
		   graphtab, comptab, output, units)

int	pcode[ARB]	# i: pseudocode used by calculator
int	maxcode		# i: maximum length of pseudocode
pointer	getvar		# i: pointer to subroutine to fetch variable values
int     nwave           # i: length of wavelength and output arrays
real    wave[ARB]       # i: wavelengths at which output is computed
char    graphtab[ARB]   # i: graph table name
char    comptab[ARB]    # i: component lookup table name
real	output[ARB]	# o: result of calculation
int	units		# o: power of FLAM (spectrum = 1, throughput = 0)
#--
int	typestack[MAXSTACK], unitstack[MAXSTACK]
pointer	valstack[MAXSTACK]

int	hi, top, bot, istr, icode, status
int	narg, ivar, ftype, ic, temp
pointer	sp, strbuf, str
real	rval
	
string	name	 "syncalc"
string	badexpr  "Synthetic photometry expression too complex"
string	badcode  "Illegal code in expression evaluator"
string	nullvar  "A variable is not legal in this expression"
string	badstat  "Error while evaluating synthetic photometry expression"
string	badfile  "Cannot read table"
string	nullexpr "Synthetic photometry expression is blank"
string  badunit  "Illegal math operation in synthetic photometry expression"

int	ctor()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (strbuf, maxcode, TY_CHAR)

	# Hi is highest array allocated on stack
	# Top is current top of stack

	hi = 0
	top = 0

	istr = 0
	icode = 1

	# Perform each operation in calculation until done

	while (pcode[icode] != Y_DONE) {

	    # Make room for results of current operation

	    top = top + 1
	    if (top > MAXSTACK)
		call synphoterr (badexpr, name)

	    if (top > hi) {
		hi = hi + 1
		call salloc (valstack[hi], nwave, TY_REAL)

	    } else if (typestack[top] == TY_STRING) {
		call salloc (valstack[top], nwave, TY_REAL)
	    }

	    switch (pcode[icode]) {
	    case Y_STR:
		# String stored in code.

		narg = 0
		icode = icode + 1
		call calcstr (pcode, icode, strbuf, istr, str)

		if (hi < MAXSTACK) {
		    hi = hi + 1
		    valstack[hi] = valstack[top]
		    typestack[hi] = typestack[top]
		}
		
		status = OK
		valstack[top] = str
		unitstack[top] = 0
		typestack[top] = TY_STRING

	    case Y_FILE:
		# Filname stored in code.

		narg = 0
		icode = icode + 1
		call calcstr (pcode, icode, strbuf, istr, str)

		typestack[top] = TY_VECTOR
		call calcread (Memc[str], nwave, wave, valstack[top], 
			       unitstack[top], status)

		if (status == ERR) {
		    if (hi < MAXSTACK) {
			hi = hi + 1
			valstack[hi] = valstack[top]
			typestack[hi] = typestack[top]
		    }
		    
		    status = OK
		    valstack[top] = str
		    unitstack[top] = 0
		    typestack[top] = TY_STRING
		}

	    case Y_NUM:
		# Number is stored in code as string. 
		# Convert and copy to stack

		ic = 1
		narg = 0
		icode = icode + 1
		call calcstr (pcode, icode, strbuf, istr, str)
		if (ctor (Memc[str], ic, rval) == 0)
		    status = ERR

		status = OK
		unitstack[top] = 0
		typestack[top] = TY_NUMBER
		call amovkr (rval, Memr[valstack[top]], nwave)

	    case Y_VAR:
		# Call user supplied routine to read variable value

		if (getvar == NULL)
		    call synphoterr (nullvar, name)

		narg = 0
		icode = icode + 1
		ivar = pcode[icode]
		call zcall2 (getvar, ivar, rval)
		
		status = OK
		unitstack[top] = 0
		typestack[top] = TY_NUMBER
		call amovkr (rval, Memr[valstack[top]], nwave)

	    case Y_FUNC:
		# Function type and number of arguments stored in code.

		icode = icode + 1
		ftype = pcode[icode]
		icode = icode + 1
		narg = pcode[icode]
		bot = top - narg
		typestack[top] = TY_VECTOR

		call calcfunc (ftype, nwave, wave, graphtab, comptab, narg, 
			       valstack[bot], typestack[bot], unitstack[bot], 
			       valstack[top], unitstack[top], status)

	    case Y_ADD:
		# Addition

		narg = 2
		typestack[top] = TY_VECTOR

		call calcadd (nwave, valstack[top-2], typestack[top-2],
			      unitstack[top-2], valstack[top], 
			      unitstack[top], status)

	    case Y_SUB:
		# Subtraction

		narg = 2
		typestack[top] = TY_VECTOR

		call calcsub (nwave, valstack[top-2], typestack[top-2],
			      unitstack[top-2], valstack[top], 
			      unitstack[top], status)

	    case Y_MUL:
		# Multiplication

		narg = 2
		typestack[top] = TY_VECTOR

		call calcmul (nwave, valstack[top-2], typestack[top-2],
			      unitstack[top-2], valstack[top], 
			      unitstack[top], status)

	    case Y_DIV:
		# Division

		narg = 2
		typestack[top] = TY_VECTOR

		call calcdiv (nwave, valstack[top-2], typestack[top-2],
			      unitstack[top-2], valstack[top], 
			      unitstack[top], status)

	    case Y_NEG:
		# Arithmetic negation

		if (typestack[top-1] == TY_STRING) {
		    status = ERR
		    call synphotwarn (badfile, Memc[valstack[top-1]])

		} else {
		    status = OK

		    narg = 1
		    typestack[top] = typestack[top-1]
		    unitstack[top] = unitstack[top-1]

		    call anegr (Memr[valstack[top-1]], 
				Memr[valstack[top]], nwave)
		}

	    default:
		call synphoterr (badcode, name)
	    }

	    # Halt if operation returns error

	    if (status == ERR)
		call synphoterr (badstat, name)

	    # Remove arguments from stack by swapping 
	    # bottom argument with result of operation

	    if (narg > 0) {
		bot = top - narg
		unitstack[bot] = unitstack[top]

		temp = valstack[bot]
		valstack[bot] = valstack[top]
		valstack[top] = temp
		
		temp = typestack[bot]
		typestack[bot] = typestack[top]
		typestack[top] = temp

		top = bot
	    }

	    icode = icode + 1
	}

	# Check to make sure there is some output

	if (top <= 0)
	    call synphoterr (nullexpr, name)

	# Check type and units of results. String type means 
	# a file that could not be read.

	if (typestack[1] == TY_STRING)
	    call synphoterr (badfile, Memc[valstack[1]])

	if (unitstack[1] < 0 || unitstack[1] > 1)
	    call synphoterr (badunit, name)

	# Copy results to output variables

	units = unitstack[1]
	call amovr (Memr[valstack[1]], output, nwave)

	call sfree (sp)

end

# CALCSTR -- Copy string from code into string buffer

procedure calcstr (pcode, icode, strbuf, istr, str)

int	pcode[ARB]	# i: pseudocode array
int	icode		# u: pointer to current code
pointer	strbuf		# i: address of string buffer
int	istr		# u: next free location in buffer
pointer	str		# o: pointer to string's location in buffer
#--

begin
	str = strbuf + istr

	while (pcode[icode] != EOS) {
	    Memc[strbuf+istr] = pcode[icode]
	    icode = icode + 1
	    istr = istr + 1
	}

	Memc[strbuf+istr] = EOS
	istr = istr + 1

end

# CALCREAD -- Read a spectrum or throughput for further processing

procedure calcread (filename, nwave, wave, value, unit, status)

char	filename[ARB]	# i: name of file containing output array
int     nwave           # i: length of wavelength and output arrays
real    wave[ARB]       # i: wavelengths at which output array is computed
pointer value		# o: pointer to output, array
int	unit		# o: units of output array
int	status		# o: error status of read operation
#--
pointer	sp, tp, wv, tabname, colname

string	wavecol  "WAVELENGTH"

int	tbcigi(), errcode()
pointer	opnsyntab()
errchk	rdband, rdspec

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (tabname, SZ_FNAME, TY_CHAR)
	call salloc (colname, SZ_FNAME, TY_CHAR)

	call breakcomp (filename, Memc[tabname], Memc[colname], SZ_FNAME)
	tp = opnsyntab (Memc[tabname])

	# If the first column is type string, what we really have is a list
	# of filenames and not a spectrum or throughput table. Return an
	# error, and postpone processing list of filenames until later.

	call syncolptr (tp, wavecol, 1, wv)
	if (tbcigi (wv, TBL_COL_DATATYPE) < 0) {
	    unit = 0
	    status = ERR

	} else {
	    # Try to open the table first as a spectrum
	    # and if that fails, as a throughput table

	    iferr {
		call rdspec (filename, nwave, wave, Memr[value])
		unit = 1
	    } then {
		# Check error code to see if error was caused by
		# column name not present in table. If not, reassert
		# error. If so, try to open as throughput table

		if (errcode () != COL_ERR) {
		    call erract (EA_ERROR)
		} else {
		    call rdband (filename, nwave, wave, Memr[value])
		    unit = 0
		}
	    }

	    status = OK
	}

	call sfree (sp)
end

# CALCFUNC -- Calculate result of function call

procedure calcfunc (ftype, nwave, wave, graphtab, comptab, narg, argval, 
		    argtype, argunit, funval, fununit, status)

int	ftype		# i: function code
int	nwave		# i: length of wavelength and ouput arrays
real    wave[ARB]       # i: wavelengths at which output is computed
char    graphtab[ARB]   # i: graph table name
char    comptab[ARB]    # i: component lookup table name
int	narg		# i: number of function arguments
pointer	argval[ARB]	# i: function arguments
int	argtype[ARB]	# i: argument types (string or vector)
int	argunit[ARB]	# i: argument units
pointer	funval		# o: output values of function
int	fununit		# o: units of output
int	status		# o: error status
#--
int	checknum[MAX_FUNC], checktype[4,MAX_FUNC], checkunit[4,MAX_FUNC]
int	i, iarg, jarg
pointer	sp, temp, func
real	poly[10]

data	checknum  / -1,  1,  2, -1,  1,  2,  4,  2,  2,  2, 
		    -1,  2,  3, -3,  4,  1,  1, -2,  2,  2  /

data	(checktype[i,1],i = 1,4)  / TY_STRING, 0, 0, 0 /
data	(checktype[i,2],i = 1,4)  / TY_NUMBER, 0, 0, 0 /
data	(checktype[i,3],i = 1,4)  / TY_NUMBER, TY_NUMBER, 0, 0 /
data	(checktype[i,4],i = 1,4)  / TY_STRING, 0, 0, 0 /
data	(checktype[i,5],i = 1,4)  / TY_NUMBER, 0, 0, 0 /
data	(checktype[i,6],i = 1,4)  / TY_NUMBER, TY_STRING, 0, 0 /
data	(checktype[i,7],i = 1,4)  / TY_NUMBER, TY_NUMBER, TY_NUMBER, TY_STRING/
data	(checktype[i,8],i = 1,4)  / TY_NUMBER, TY_NUMBER, 0, 0 /
data	(checktype[i,9],i = 1,4)  / TY_STRING, TY_NUMBER, 0, 0 /
data	(checktype[i,10],i = 1,4) / TY_NUMBER, TY_NUMBER, 0, 0 /
data	(checktype[i,11],i = 1,4) / TY_STRING, 0, 0, 0 /
data	(checktype[i,12],i = 1,4) / TY_NUMBER, TY_NUMBER, 0, 0 /
data	(checktype[i,13],i = 1,4) / TY_NUMBER, TY_NUMBER, TY_STRING, 0 /
data	(checktype[i,14],i = 1,4) / TY_NUMBER, TY_NUMBER, TY_NUMBER, 0 /
data	(checktype[i,15],i = 1,4) / TY_VECTOR, TY_VECTOR, TY_NUMBER, TY_STRING/
data	(checktype[i,16],i = 1,4) / TY_STRING, 0, 0, 0 /
data	(checktype[i,17],i = 1,4) / TY_STRING, 0, 0, 0 /
data	(checktype[i,18],i = 1,4) / TY_VECTOR, TY_NUMBER, 0, 0 /
data	(checktype[i,19],i = 1,4) / TY_NUMBER, TY_STRING, 0, 0 /
data	(checktype[i,20],i = 1,4) / TY_VECTOR, TY_NUMBER, 0, 0 /

data	(checkunit[i,1],i = 1,4)  / 0, 0, 0, 0 /
data	(checkunit[i,2],i = 1,4)  / 0, 0, 0, 0 /
data	(checkunit[i,3],i = 1,4)  / 0, 0, 0, 0 /
data	(checkunit[i,4],i = 1,4)  / 0, 0, 0, 0 /
data	(checkunit[i,5],i = 1,4)  / 0, 0, 0, 0 /
data	(checkunit[i,6],i = 1,4)  / 0, 0, 0, 0 /
data	(checkunit[i,7],i = 1,4)  / 0, 0, 0, 0 /
data	(checkunit[i,8],i = 1,4)  / 0, 0, 0, 0 /
data	(checkunit[i,9],i = 1,4)  / 0, 0, 0, 0 /
data	(checkunit[i,10],i = 1,4) / 0, 0, 0, 0 /
data	(checkunit[i,11],i = 1,4) / 0, 0, 0, 0 /
data	(checkunit[i,12],i = 1,4) / 0, 0, 0, 0 /
data	(checkunit[i,13],i = 1,4) / 0, 0, 0, 0 /
data	(checkunit[i,14],i = 1,4) / 0, 0, 0, 0 /
data	(checkunit[i,15],i = 1,4) / 1, 0, 0, 0 /
data	(checkunit[i,16],i = 1,4) / 0, 0, 0, 0 /
data	(checkunit[i,17],i = 1,4) / 0, 0, 0, 0 /
data	(checkunit[i,18],i = 1,4) / 0, 0, 0, 0 /
data	(checkunit[i,19],i = 1,4) / 0, 0, 0, 0 /
data	(checkunit[i,20],i = 1,4) / 1, 0, 0, 0 /

string	funcstr  FUNCNAMES
string	catfile  CAT_TEMPLATE
string	badnum   "Incorrect number of arguments for function"
string	badfile  "Cannot access file from function"
string	badtype  "Incorrect argument type for function"
string	badunit  "Incorrect argument units for function"

int	word_find()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (temp, SZ_LINE, TY_CHAR)
	call salloc (func, SZ_FNAME, TY_CHAR)

	# Get function name for error messages

	if (word_find (ftype, funcstr, Memc[func], SZ_FNAME) == 0)
	    call strcpy ("unknown", Memc[func], SZ_FNAME)

	status = OK

	# Check number of arguments

	if (checknum[ftype] > 0 && checknum[ftype] != narg) {
	    call synphotwarn (badnum, Memc[func])
	    status = ERR
	}

	# Check type and units of arguments

	jarg = min (narg, abs (checknum[ftype]))
	do iarg = 1, jarg {
	    if (checktype[iarg,ftype] != argtype[iarg]) {
		if (checktype[iarg,ftype] == TY_VECTOR &&
		    argtype[iarg] == TY_STRING) {
		    call synphotwarn (badfile, Memc[func])

		} else {
		    call synphotwarn (badtype, Memc[func])
		}
		status = ERR

	    } else if (checkunit[iarg,ftype] != argunit[iarg]) {
		call synphotwarn (badunit, Memc[func])
		status = ERR
	    }
	}

	# Check types and units of function with an indeterminant
	# number of arguments. Here the type and units of the last
	# fixed argument is the same as the type and units of the 
	# remaining arguments.

	if (checknum[ftype] < 0) {
	    jarg = abs (checknum[ftype])
	    do iarg = jarg+1, narg {
		if (checktype[jarg,ftype] != argtype[iarg]) {
		    call synphotwarn (badtype, Memc[func])
		    status = ERR

		} else if (checkunit[jarg,ftype] != argunit[iarg]) {
		    call synphotwarn (badunit, Memc[func])
		    status = ERR
		}
	    }
	}

	# Call the appropriate function to calculate the results

	if (status == ERR) {
	    funval = NULL
	    fununit = 0

	} else {
	    # Concatenate string arguments

	    if (ftype == FN_BAND || ftype == FN_CAT || ftype == FN_ICAT)
		call mkmodestr (argval, narg, Memc[temp], SZ_LINE)

	    switch (ftype) {
	    case FN_BAND:	# instrument bandpass 
		fununit = 0
		call bandfunc (Memc[temp], graphtab, comptab, 
			       nwave, wave, Memr[funval])
		
	    case FN_BB:	# black body spectrum
		fununit = 1
		call bbfunc (Memr[argval[1]], nwave, wave, Memr[funval])
		call amulkr (Memr[funval], RENORM, Memr[funval], nwave)
		
	    case FN_BOX:	# box function
		fununit = 0
		call boxfunc (Memr[argval[1]], Memr[argval[2]], 
			      nwave, wave, Memr[funval])
		
	    case FN_CAT:
		fununit = 1
		call catfunc (Memc[temp], catfile, nwave, wave, Memr[funval])

	    case FN_EBMV:	# galactic extinction function
		fununit = 0
		call ebmvfunc (Memr[argval[1]], nwave, wave, Memr[funval])
		
	    case FN_EBMVX:	# extended extinction function
		fununit = 0
		call ebmvxfunc (Memr[argval[1]], Memc[argval[2]],
				nwave, wave, Memr[funval])
		
	    case FN_EM:		# emission line function
		fununit = 1
		call emfunc (Memr[argval[1]], Memr[argval[2]], 
			     Memr[argval[3]], Memc[argval[4]],
			     nwave, wave, Memr[funval])

	    case FN_GAUSS:	# gaussian function
		fununit = 0
		call gaussfunc (Memr[argval[1]], Memr[argval[2]], 
				nwave, wave, Memr[funval])
		
	    case FN_GRID:	# interpolation on grid of spectra
		fununit = 1
		call gridfunc (Memc[argval[1]], Memr[argval[2]],
			       nwave, wave, Memr[funval])
		
	    case FN_HI:		# hydrogen emission spectrum
		fununit = 1
		call hifunc (Memr[argval[1]], Memr[argval[2]], 
			     nwave, wave, Memr[funval])
		call amulkr (Memr[funval], RENORM, Memr[funval], nwave)
		
	    case FN_ICAT:
		fununit = 1
		call icatfunc (Memc[temp], catfile, nwave, wave, Memr[funval])

	    case FN_LGAUSS:	# log gauss function
		fununit = 0
		call lgaussfunc (Memr[argval[1]], Memr[argval[2]], 
				 nwave, wave, Memr[funval])
		
	    case FN_PL:		# power law spectrum
		fununit = 1
		call plfunc (Memr[argval[1]], Memr[argval[2]], Memc[argval[3]],
			     nwave, wave, Memr[funval])

	    case FN_POLY:	# legendre polynomial
		jarg = min (narg, 10+2)
		do iarg = 3, jarg
		    poly[iarg-2] = Memr[argval[iarg]]
		
		fununit = 0
		call polyfunc (Memr[argval[1]], Memr[argval[2]], jarg-2, poly,
			       nwave, wave, Memr[funval])

	    case FN_RN:		# renormalization
		fununit = 1
		call renormfunc (Memr[argval[1]], Memr[argval[2]], 
				 Memr[argval[3]], Memc[argval[4]],
				 nwave, wave, Memr[funval])

	    case FN_SPEC:	# read a spectrum
		fununit = 1
		call rdspec (Memc[argval[1]], nwave, wave, Memr[funval])

	    case FN_THRU:	# read a throughput table
		fununit = 0
		call rdband (Memc[argval[1]], nwave, wave, Memr[funval])

	    case FN_TILT:	# modify passband by legendre polynomial
		jarg = min (narg, 10+1)
		do iarg = 2, jarg
		    poly[iarg-1] = Memr[argval[iarg]]
		
		fununit = 0
		call tiltfunc (Memr[argval[1]], jarg-1, poly, 
			       nwave, wave, Memr[funval])

	    case FN_UNIT:	# constant spectrum
		fununit = 1
		call unitfunc (Memr[argval[1]], Memc[argval[2]],
			       nwave, wave, Memr[funval])

	    case FN_Z:		# redshift spectrum
		fununit = 1
		call zfunc (Memr[argval[1]], Memr[argval[2]],
			    nwave, wave, Memr[funval])
	    }
	}

	call sfree (sp)
end

# CALCADD -- Add two vectors to produce output vector

procedure calcadd (nwave, argval, argtype, argunit, outval, outunit, status)

int	nwave		# i: length of wavelength and ouput arrays
pointer	argval[2]	# i: operation arguments
int	argtype[2]	# i: argument types (string or vector)
int	argunit[2]	# i: argument units
pointer	outval		# o: output values of function
int	outunit		# o: units of output
int	status		# o: error status
#--
int	iarg

string	badfile  "Cannot read table"
string	badunit  "Units must match to perform addition"

begin
	status = OK
	do iarg = 1, 2 {
	    if (argtype[iarg] == TY_STRING) {
		call synphotwarn (badfile, Memc[argval[iarg]])
		status = ERR
	    }
	}

	outunit = argunit[1]
	if (argunit[1] != argunit[2]) {
	    call synphotwarn (badunit, "syncalc")
	    status = ERR
	}

	if (status != ERR)
	    call aaddr (Memr[argval[1]], Memr[argval[2]], Memr[outval], nwave)

end

# CALCSUB -- Subtract two vectors to produce output vector

procedure calcsub (nwave, argval, argtype, argunit, outval, outunit, status)

int	nwave		# i: length of wavelength and ouput arrays
pointer	argval[2]	# i: operation arguments
int	argtype[2]	# i: argument types (string or vector)
int	argunit[2]	# i: argument units
pointer	outval		# o: output values of function
int	outunit		# o: units of output
int	status		# o: error status
#--
int	iarg

string	badfile  "Cannot read table"
string	badunit  "Units must match to perform subtraction"

begin
	status = OK
	do iarg = 1, 2 {
	    if (argtype[iarg] == TY_STRING) {
		call synphotwarn (badfile, Memc[argval[iarg]])
		status = ERR
	    }
	}

	outunit = argunit[1]
	if (argunit[1] != argunit[2]) {
	    call synphotwarn (badunit, "syncalc")
	    status = ERR
	}

	if (status != ERR)
	    call asubr (Memr[argval[1]], Memr[argval[2]], Memr[outval], nwave)

end

# CALCMUL -- Multiply two vectors to produce output vector

procedure calcmul (nwave, argval, argtype, argunit, outval, outunit, status)

int	nwave		# i: length of wavelength and ouput arrays
pointer	argval[2]	# i: operation arguments
int	argtype[2]	# i: argument types (string or vector)
int	argunit[2]	# i: argument units
pointer	outval		# o: output values of function
int	outunit		# o: units of output
int	status		# o: error status
#--
int	iarg

string	badfile  "Cannot read table"

begin
	status = OK
	do iarg = 1, 2 {
	    if (argtype[iarg] == TY_STRING) {
		call synphotwarn (badfile, Memc[argval[iarg]])
		status = ERR
	    }
	}

	if (status != ERR) {
	    call amulr (Memr[argval[1]], Memr[argval[2]], Memr[outval], nwave)
	    outunit = argunit[1] + argunit[2]
	}

end

# CALCDIV -- Divide two vectors to produce output vector

procedure calcdiv (nwave, argval, argtype, argunit, outval, outunit, status)

int	nwave		# i: length of wavelength and ouput arrays
pointer	argval[2]	# i: operation arguments
int	argtype[2]	# i: argument types (string or vector)
int	argunit[2]	# i: argument units
pointer	outval		# o: output values of function
int	outunit		# o: units of output
int	status		# o: error status
#--
int	iarg
real    calczero()
extern	calczero

string	badfile  "Cannot read table"

begin
	status = OK
	do iarg = 1, 2 {
	    if (argtype[iarg] == TY_STRING) {
		call synphotwarn (badfile, Memc[argval[iarg]])
		status = ERR
	    }
	}

	if (status != ERR) {
	    call advzr (Memr[argval[1]], Memr[argval[2]], 
			Memr[outval], nwave, calczero)
	    outunit = argunit[1] - argunit[2]
	}

end

# CALCZERO -- Zero replacement for division procedure

real procedure calczero (value)

real	value
#--

begin
	return (0.0)
end

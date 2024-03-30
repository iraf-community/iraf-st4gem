include	<mach.h>
include	<tbset.h>
include	<synphot.h>
include "syncalc.h"
include	"libsynphot.h"

# History:
# H.Bushouse 1 Nov. 2002: syw_addwave modified to initialize variable iw=0.
#			  syw_copy modified to change use of amovr to amovi
#			  for copying SW_KEEPER array.
#			  syw_func modified FN_Z case so that only the max
#			  wave limit gets increased for z>0 and only min wave
#			  limit gets decreased for z<0.

define	MINLENGTH	100
define	MIN_DELTA	(80.0*EPSILONR)

define	MAXSTACK	64
define	FUNCNAME	"synwave"

define	LEN_RGSTRUCT	5

define	SW_STRING	1
define	SW_NUMBER	2
define	SW_LIST		3

define	SW_TYPE		Memi[$1]
define	SW_LENGTH	Memi[$1+1]
define	SW_LOPTR	Memi[$1+2]
define	SW_HIPTR	Memi[$1+3]
define	SW_KPTR		Memi[$1+4]

define	SW_STR		Memc[SW_LOPTR($1)]		
define	SW_VAL		Memr[SW_LOPTR($1)]
define	SW_LORANGE	Memr[SW_LOPTR($1)+($2)]
define	SW_HIRANGE	Memr[SW_HIPTR($1)+($2)]
define	SW_KEEPER	Memi[SW_KPTR($1)+($2)]

define	SW_ERROR	(SW_LENGTH($1) == 0)

# SYNWAVE -- Calculate wavelength range appropriate for expression

procedure synwave (intersect, pcode, ncode, maxcode, graphtab, comptab, 
		   wave, nwave)

int	intersect		# i: intersect=YES union=NO
int	pcode[maxcode,ncode]	# i: Expression pseudocode arrays
int	ncode			# i: Number of pseudocode arrays
int	maxcode			# i: dimension of pseudocode array
char	graphtab[ARB]		# i: Instrument graph table
char	comptab[ARB]		# i: Component name table
pointer	wave			# o: wavelength set (must be freed by caller)
int	nwave			# o: length of wavelength set
#--
int	icode, status
pointer	rg1, rg2, rg3
real	minwave, maxwave

string	badvega  "Error in Vega spectrum"

int	getnwave()
pointer	syw_calc(), syw_create(), syw_union(), syw_intersect()

begin

	nwave = LENWAVE
	do icode = 1, ncode {
	    nwave = max (nwave, getnwave (pcode[1,icode], maxcode))
	    rg2 = syw_calc (pcode[1,icode], maxcode, graphtab, comptab)

	    if (icode == 1) {
		rg1 = rg2

	    } else {
		if (intersect == NO) {
		    rg3 = syw_union (rg1, rg2)
		} else {
		    rg3 = syw_intersect (rg1, rg2)
		}

		call syw_free (rg1)
		call syw_free (rg2)
		rg1 = rg3
	    }
	}

	# Replace INDEF values with defaults

	call syw_bound (rg1, DEFMIN, DEFMAX)

	# Take intersection with vega's spectrum

	call tabrange (VEGA, minwave, maxwave, status)
	if (status == ERR)
	    call synphoterr (badvega, VEGA)

	rg2 = syw_create (minwave, maxwave)
	rg3 = syw_intersect (rg1, rg2)

	call syw_free (rg1)
	call syw_free (rg2)

	call syw_wave (rg3, wave, nwave)
	call syw_free (rg3)
end

# SYW_ADDWAVE -- Append a new wavelength set to the output set

procedure syw_addwave (logspace, rlo, rhi, nwave, maxlen, iwave, wave)
		       

bool	logspace	# i: Use logarithmic spacing for wavelengths?
real	rlo		# i: low endpoint of range
real	rhi		# i: high endpoint of range
int	nwave		# i: number of points in wavelength set
int	maxlen		# i: allocated length of wave set
int	iwave		# u: current length of wavelength set
pointer	wave		# i: wavelength set
#--
int	iw, mwave
pointer	sp, wpart

begin
	# Check for boundary case

	if (nwave <= 0)
	    return

	# Allocate temporary array for wavelengths

	call smark (sp)
	call salloc (wpart, nwave, TY_REAL)

	# Create the wavelength set

	call waveset (logspace, rlo, rhi, nwave, Memr[wpart])

	# Delete wavelengths that overlap previous wavelengths
        iw = 0
	if (iwave > 0) {
	    for (iw = 0; iw < nwave; iw = iw + 1) {
		if (Memr[wpart+iw] > Memr[wave+iwave-1])
		    break
	    }
	}

	# Compute number of wavelengths to copy

	mwave = nwave - iw
	if (iwave + mwave > maxlen)
	    mwave = maxlen - iwave

	# Copy wavelengths into output array and update index

	if (mwave > 0) {
	    call amovr (Memr[wpart+iw], Memr[wave+iwave], mwave)
	    iwave = iwave + mwave
	}

	call sfree (sp)
end

# SYW_BOUND -- Replace indef values with substitute values

procedure syw_bound (rg, minwave, maxwave)

pointer	rg		# u: the range descriptor
real	minwave		# i: replacement for low indef values
real	maxwave		# i: replacement for high indef values
#--
int	ir

begin
	do ir = 0, SW_LENGTH(rg)-1 {
	    if (IS_INDEFR(SW_LORANGE(rg,ir)))
		SW_LORANGE(rg,ir) = minwave

	    if (IS_INDEFR(SW_HIRANGE(rg,ir)))
		SW_HIRANGE(rg,ir) = maxwave
	}

end

# SYW_CALC -- Calculate  a wavelength range set for a single expression

pointer procedure syw_calc (pcode, maxcode, graphtab, comptab)

int	pcode[ARB]	# i: pseudocode used by calculator
int	maxcode		# i: maximum length of pseudocode
char    graphtab[ARB]   # i: graph table name
char    comptab[ARB]    # i: component lookup table name
#--
pointer	stack[MAXSTACK]

int	top, bot, icode, ftype, iarg, narg
pointer	sp, strbuf, rbuf, temp

string	badexpr  "Synthetic photometry expression too complex"
string	badcode  "Illegal code in expression evaluator"
string	nullvar  "A variable is not legal in this expression"
string	badstat  "Error while evaluating synthetic photometry expression"
string	badfile  "Cannot read table"
string	nullexpr "Synthetic photometry expression is blank"

pointer	syw_make(), syw_func(), syw_math(), syw_default()
pointer	syw_negate(), syw_table()
real	syw_real()
int	strlen()

begin
	# Allocate temporary buffers

	call smark (sp)
	call salloc (strbuf, SZ_FNAME, TY_CHAR)
	call salloc (rbuf, 1, TY_REAL)

	# Loop over instructions in pseudocode

	top = 0
	for (icode = 1; pcode[icode] != Y_DONE; icode = icode + 1) {
	    # Make room for results of current operation

	    top = top + 1
	    if (top > MAXSTACK)
		call synphoterr (badexpr, FUNCNAME)

	    switch (pcode[icode]) {
	    case Y_STR:
		# Push a string on the stack

		narg = 0
		icode = icode + 1
		call syw_str (pcode, icode, Memc[strbuf], SZ_FNAME)
		stack[top] = syw_make (strbuf, SZ_FNAME, SW_STRING) 

	    case Y_FILE:
		# Read the wavelength set of a spectrum

		narg = 0
		icode = icode + 1
		call syw_str (pcode, icode, Memc[strbuf], SZ_FNAME)
		stack[top] = syw_table (Memc[strbuf])

		# Couldn't open file name must actually be string

		if (SW_ERROR(stack[top])) {
		    call syw_free (stack[top])
		    stack[top] = syw_make (strbuf, SZ_FNAME, SW_STRING) 
		}
		
	    case Y_NUM:
		# Push a number on the stack

		narg = 0
		icode = icode + 1
		call syw_str (pcode, icode, Memc[strbuf], SZ_FNAME)
		Memr[rbuf] = syw_real (Memc[strbuf])
		stack[top] = syw_make (rbuf, 1, SW_NUMBER)
		
	    case Y_VAR:
		# Variable zero

		narg = 0
		icode = icode + 1
		Memr[rbuf] = INDEFR
		stack[top] = syw_make (rbuf, 1, SW_NUMBER)

	    case Y_FUNC:
		# Evaluate a built in function

		icode = icode + 1
		ftype = pcode[icode]
		icode = icode + 1
		narg = pcode[icode]
		bot = top - narg

		stack[top] = syw_func (ftype, graphtab ,comptab, 
				       narg, stack[bot])

	    case Y_ADD, Y_SUB, Y_MUL, Y_DIV:
		# Take union of two ranges

		narg = 2
		stack[top] = syw_math (pcode[icode], stack[top-2], 
				       stack[top-1])
 
	    case Y_NEG:
		# Negation

		narg = 1
		stack[top] = syw_negate (stack[top-1])

	    default:
		# Illegal instruction

		call synphoterr (badcode, FUNCNAME)
	    }

	    if (SW_ERROR(stack[top])) {
		call synphoterr (badstat, FUNCNAME)
	    }

	    # Pop arguments off top of stack and release memory

	    if (narg > 0) {
		bot = top - narg

		temp = stack[bot]
		stack[bot] = stack[top]
		stack[top] = temp

		do iarg = 1, narg
		    call syw_free (stack[bot+iarg])

		top = bot
	    }
	}

	# Check to make sure there is some output

	if (top <= 0)
	    call synphoterr (nullexpr, FUNCNAME)

	# Check type and length of results. String type means 
	# a file that could not be read.

	if (SW_TYPE(stack[1]) == SW_STRING) {
	    call strcpy (SW_STR(stack[1]),Memc[strbuf],strlen(SW_STR(stack[1])))
	    call syw_free (stack[1])
	    call synphoterr (badfile, Memc[strbuf])
	}

	if (SW_TYPE(stack[1]) == SW_NUMBER) {
	    call syw_free (stack[1])
	    stack[1] = syw_default ()
	}

	return (stack[1])
end

# SYW_COMPARE -- Compare the ordering of two ranges

int procedure syw_compare (rg, ix, jx)

pointer	rg		# i: wavelength set descriptor
int	ix		# i: first index into array
int	jx		# i: second index into array
#--
int	order
real	r1, r2, q1, q2

begin
	# Translate indices into range values

	if (ix > 0) {
	    r1 = SW_LORANGE(rg,ix-1)
	    q1 = SW_HIRANGE(rg,ix-1)
	} else {
	    r1 = SW_HIRANGE(rg,-ix-1)
	    q1 = SW_LORANGE(rg,-ix-1)
	}

	if (jx > 0) {
	    r2 = SW_LORANGE(rg,jx-1)
	    q2 = SW_HIRANGE(rg,jx-1)
	} else {
	    r2 = SW_HIRANGE(rg,-jx-1)
	    q2 = SW_LORANGE(rg,-jx-1)
	}

	if (r1 < r2) {
	    order = -1
	} else if (r1 > r2) {
	    order = 1
	} else if (q1 < q2) {
	    order = 1
	} else if (q1 > q2) {
	    order = -1
	} else {
	    order = 0
	}

	return (order)
end

# SYW_COPY -- Copy a wavelength set descriptor

pointer procedure syw_copy (rg1)

pointer	rg1		# i: wavelength set descriptor
#--
int	nr
pointer	rg2

pointer	syw_new()

begin
	nr = SW_LENGTH(rg1)
	rg2 = syw_new (nr, SW_TYPE(rg1))

	if (SW_LOPTR(rg1) != NULL)
	    call amovr (SW_LORANGE(rg1,0), SW_LORANGE(rg2,0), nr)

	if (SW_HIPTR(rg1) != NULL)
	    call amovr (SW_HIRANGE(rg1,0), SW_HIRANGE(rg2,0), nr)

	if (SW_KPTR(rg1) != NULL)
	    call amovi (SW_KEEPER(rg1,0), SW_KEEPER(rg2,0), nr)

	return (rg2)
end

# SYW_CREATE -- Create a wavelength range with specified endpoints

pointer procedure syw_create (minwave, maxwave)

real	minwave		# i: start of wavelength range
real	maxwave		# i: end of wavelength range
#--
pointer	rg

pointer	syw_new()

begin
	# Allocate and initialize structure

	rg = syw_new (1, SW_LIST)

	SW_LORANGE(rg,0) = minwave
	SW_HIRANGE(rg,0) = maxwave
	SW_KEEPER(rg,0) = NO

	return (rg)
end

# SYW_DEFAULT -- Create a default wavelength range

pointer	procedure syw_default ()

#--
pointer	syw_create()

begin
	# INDEF values indicate the endpoints 
	# of the range are undefined

	return syw_create (INDEFR, INDEFR)
end

# SYW_ERROR -- Create an error range

pointer procedure syw_error ()

#--
pointer	syw_new()

begin
	# Create an empty structure to signal error

	return (syw_new (0, SW_LIST))
end

# SYW_FREE -- Free a wave descriptor

procedure syw_free (rg)

pointer	rg		# i: wave descriptor
#--

begin
	if (SW_LOPTR(rg) != NULL) {
	    if (SW_TYPE(rg) == SW_STRING) {
		call mfree (SW_LOPTR(rg), TY_CHAR)
	    } else {
		call mfree (SW_LOPTR(rg), TY_REAL)
	    }
	}

	if (SW_HIPTR(rg) != NULL)
	    call mfree (SW_HIPTR(rg), TY_REAL)

	if (SW_KPTR(rg) != NULL)
	    call mfree (SW_KPTR(rg), TY_REAL)

	call mfree (rg, TY_INT)
end

# SYW_FUNC -- Compute the wavelength set for a function

pointer procedure syw_func (ftype, graphtab, comptab, narg, stack)

int	ftype		# i: pseudocode representation of function type
char    graphtab[ARB]   # i: graph table name
char    comptab[ARB]    # i: component lookup table name
int	narg		# i: number of function arguments
int	stack[ARB]	# i: function arguments as descriptors
#--
int	nfile, interp, ir
pointer	sp, mode, filelist, rg
real	minwave, maxwave, factor

string	catfile  CAT_TEMPLATE
string	badfile  "Cannot read table"
string	badtype  "Incorrect argument type for function"

int	btoi()
pointer	syw_table(), syw_copy(), syw_list(), syw_default(), syw_create()
pointer	syw_union()

begin
	call smark (sp)
	call salloc (mode, SZ_LINE, TY_CHAR)
	call salloc (filelist, (SZ_FNAME+1)*MAXLIST, TY_CHAR)

	switch (ftype) {
	case FN_BAND:
	    # Bandpass function

	    call syw_mode (stack, narg, Memc[mode], SZ_LINE)

	    call graffiles (Memc[mode], graphtab, comptab, SZ_FNAME, 
			    MAXLIST, nfile, Memc[filelist])

	    rg = syw_list (SZ_FNAME, nfile, Memc[filelist])

	case FN_BOX:
	    # Box function

	    if (SW_TYPE(stack[1]) != SW_NUMBER || 
		SW_TYPE(stack[2]) != SW_NUMBER   )
		call synphoterr (badtype, FUNCNAME)

	    if (IS_INDEFR(SW_VAL(stack[1])) || IS_INDEFR(SW_VAL(stack[2]))) {
		rg = syw_default ()

	    } else {
		minwave = SW_VAL(stack[1]) - 0.5 * SW_VAL(stack[2])
		maxwave = SW_VAL(stack[1]) + 0.5 * SW_VAL(stack[2])
		rg = syw_create (minwave, maxwave)
	    }

	case FN_CAT, FN_ICAT:
	    # Interpolate in wavelength catalog

	    interp = btoi (ftype == FN_ICAT)
	    call syw_mode (stack, narg, Memc[mode], SZ_LINE)

	    call catfiles (Memc[mode], catfile, interp, SZ_FNAME,
			   MAXLIST, nfile, Memc[filelist])

	    rg = syw_list (SZ_FNAME, nfile, Memc[filelist])

	case FN_EM, FN_GAUSS:
	    # Emission line or gaussian passband

	    if (SW_TYPE(stack[1]) != SW_NUMBER || 
		SW_TYPE(stack[2]) != SW_NUMBER)
		call synphoterr (badtype, FUNCNAME)

	    if (IS_INDEFR(SW_VAL(stack[1])) || IS_INDEFR(SW_VAL(stack[2]))) {
		rg = syw_default ()

	    } else {
		factor = 6.0 / sqrt (8.0 * log (2.0))
		minwave = SW_VAL(stack[1]) - factor * SW_VAL(stack[2])
		maxwave = SW_VAL(stack[1]) + factor * SW_VAL(stack[2])

		rg = syw_create (minwave, maxwave)
	    }

	case FN_GRID:
	    # Interpolate between spectra in file

	    if (SW_TYPE(stack[1]) != SW_STRING || 
		SW_TYPE(stack[2]) != SW_NUMBER)
		call synphoterr (badtype, FUNCNAME)

	    if (IS_INDEFR(SW_VAL(stack[2]))) {
		rg = syw_default ()

	    } else {
		call gridfiles (SW_STR(stack[1]), SW_VAL(stack[2]), 
				SZ_FNAME, MAXLIST, nfile, Memc[filelist])

		rg = syw_list (SZ_FNAME, nfile, Memc[filelist])
	    }

	case FN_RN:
	    # Renormalize spectrum
	    rg = stack[2]
	    do ir = 0, SW_LENGTH(rg)-1
		SW_KEEPER(rg,ir) = YES

	    rg = syw_union (stack[1], stack[2])

	case FN_SPEC, FN_THRU:
	    # Spectrum file

	    rg = syw_table (SW_STR(stack[1]))
	    if (SW_ERROR (rg))
		call synphoterr (badfile, SW_STR(stack[1]))

	case FN_Z:
	    # Red shift

	    if (SW_TYPE(stack[2]) != SW_NUMBER)
		call synphoterr (badtype, FUNCNAME)

	    if (IS_INDEFR(SW_VAL(stack[2]))) {
		rg = syw_default ()

	    } else {
	        rg = syw_copy (stack[1])
	        factor = 1.0 + SW_VAL(stack[2])

	        do ir = 0, SW_LENGTH(rg)-1 {
		   if (SW_KEEPER(rg,ir) == YES)
		       next

		   if (IS_INDEFR(SW_LORANGE(rg,ir)) ||
		       IS_INDEFR(SW_HIRANGE(rg,ir)))
		       next

		   if (factor < 1.0)
		       SW_LORANGE(rg,ir) = factor * SW_LORANGE(rg,ir)
		   else
		       SW_HIRANGE(rg,ir) = factor * SW_HIRANGE(rg,ir)
	        }
	    }


	default:
	    rg= syw_default ()
	}

	call sfree (sp)
	return (rg)
end

# SYW_INTERSECT -- Take the intersection of two range sets

pointer procedure syw_intersect (rg1, rg2)

pointer	rg1	# i: first descriptor
pointer	rg2	# i: second descriptor
#--
int	ir, jr, kr, nrange
pointer	rg3

string	badtype  "syw_intersect: can only intersect lists"

pointer	syw_new()

begin
	# Check descriptors to make sure they are wavelength sets

	if (SW_TYPE(rg1) != SW_LIST || SW_TYPE(rg2) != SW_LIST)
	    call synphoterr (badtype, FUNCNAME)

	# Create new wavelength descriptor

	nrange = SW_LENGTH(rg1) * SW_LENGTH(rg2)
	rg3 = syw_new (nrange, SW_LIST)
	
	# Intersect the ranges from the two descriptors

	kr = 0
	for (ir=0; ir < SW_LENGTH(rg1); ir=ir+1) {
	    if (SW_KEEPER(rg1,ir) == YES)
		next

	    for (jr=0; jr < SW_LENGTH(rg2); jr=jr+1) {
		if (SW_KEEPER(rg2,jr) == YES)
		    next

		if (IS_INDEFR(SW_LORANGE(rg1,ir)) && 
		    IS_INDEFR(SW_HIRANGE(rg1,ir))) {

		    SW_LORANGE(rg3,kr) = SW_LORANGE(rg2,jr)
		    SW_HIRANGE(rg3,kr) = SW_HIRANGE(rg2,jr)
		    SW_KEEPER(rg3,kr) = NO
		    kr = kr + 1

		} else if (IS_INDEFR(SW_LORANGE(rg2,jr)) && 
			   IS_INDEFR(SW_HIRANGE(rg2,jr))) {

		    SW_LORANGE(rg3,kr) = SW_LORANGE(rg1,ir)
		    SW_HIRANGE(rg3,kr) = SW_HIRANGE(rg1,ir)
		    SW_KEEPER(rg3,kr) = NO
		    kr = kr + 1

		} else if (SW_HIRANGE(rg1,ir) > SW_LORANGE(rg2,jr) &&
			   SW_LORANGE(rg1,ir) < SW_HIRANGE(rg2,jr)   ) {

		    SW_LORANGE(rg3,kr) = max (SW_LORANGE(rg1,ir), 
					      SW_LORANGE(rg2,jr))
		    SW_HIRANGE(rg3,kr) = min (SW_HIRANGE(rg1,ir),
					      SW_HIRANGE(rg2,jr))
		    SW_KEEPER(rg3,kr) = NO
		    kr = kr + 1
		}
	    }
	}

	# Copy the kept ranges to the output descriptor

	for (ir=0; ir < SW_LENGTH(rg1); ir=ir+1) {
	    if (SW_KEEPER(rg1,ir) == NO)
		next

	    SW_LORANGE(rg3,kr) = SW_LORANGE(rg1,ir)
	    SW_HIRANGE(rg3,kr) = SW_HIRANGE(rg1,ir)
	    SW_KEEPER(rg3,kr) = YES
	    kr = kr + 1
	}

	for (jr=0; jr < SW_LENGTH(rg2); jr=jr+1) {
	    if (SW_KEEPER(rg2,jr) == NO)
		next

	    SW_LORANGE(rg3,kr) = SW_LORANGE(rg2,jr)
	    SW_HIRANGE(rg3,kr) = SW_HIRANGE(rg2,jr)
	    SW_KEEPER(rg3,kr) = YES
	    kr = kr + 1
	}

	SW_LENGTH(rg3) = kr
	return (rg3)
end

# SYW_LIMITS -- Compute the minimum and maximum value from a range descriptor

procedure syw_limits (rg, minwave, maxwave)

pointer	rg			# i: range descriptor
real	minwave			# o: lower limit of range
real	maxwave			# o: upper limit of range
#--
int	ir

begin
	if (SW_LENGTH(rg) < 1) {
	    minwave = DEFMIN
	    maxwave = DEFMAX
	} else {
	    minwave = SW_LORANGE(rg,0)
	    maxwave = SW_HIRANGE(rg,0)
	}

	do ir = 1, SW_LENGTH(ir) - 1 {
	    minwave = min (minwave, SW_LORANGE(rg,ir))
	    maxwave = max (maxwave, SW_HIRANGE(rg,ir))
	}

end

# SYW_LIST -- Compute the wavelength range of a list of files

pointer procedure syw_list (maxfile, nfile, filelist)

int	maxfile			# i: maximum length of file name
int	nfile			# i: number of instrument components
char	filelist[maxfile,ARB]	# i: list of component file names
#--
int	ifile
pointer	rg1, rg2, rg3

errchk	tabrange
pointer	syw_error(), syw_table(), syw_intersect()

begin
	if (nfile <= 0) {
	    rg1 = syw_error ()

	} else {
	    # Loop over all spectra

	    do ifile = 1, nfile {
		rg2 = syw_table (filelist[1,ifile])

		if (SW_ERROR (rg2)) {
		    if (ifile > 1)
			call syw_free (rg1)

		    rg1 = rg2
		    break
		}

		# Take intersection of range with previous ranges

		if (ifile == 1) {
		    rg1 = rg2

		} else {
		    rg3 = syw_intersect (rg1, rg2)

		    call syw_free (rg1)
		    call syw_free (rg2)
		    rg1 = rg3
		}    
	    }
	}

	return (rg1)
end

# SYW_MAKE -- Copy an array into a range descriptor

pointer	procedure syw_make (buffer, length, type)

pointer	buffer		# i: the array to copy
int	length		# i: the length of the array
int	type		# i: the type of the array
#--
int	rg

string	badtype  "syw_make: invalid input type"

pointer	syw_new()

begin
	rg = syw_new (length, type)

	switch (type) {
	case SW_STRING:
	    call amovc (Memc[buffer], Memc[SW_LOPTR(rg)], length)

	case SW_NUMBER:
	    call amovr (Memr[buffer], Memr[SW_LOPTR(rg)], length)

	default:
	    call synphoterr (badtype, FUNCNAME)
	}

	return (rg)
end

# SYW_MATH -- Perform a binary arithmetic operation on two descriptors

pointer procedure syw_math (op, rg1, rg2)

int	op		# i: operation to perform
pointer	rg1		# i: first descriptor
pointer	rg2		# i: second descriptor
#--
pointer	rg3

string	badfile  "Cannot read table"

pointer	syw_copy(), syw_union(), syw_intersect()

begin
	# A string value indicates that a file was not readable

	if (SW_TYPE(rg1) == SW_STRING)
	    call synphoterr (badfile, SW_STR(rg1))

	if (SW_TYPE(rg2) == SW_STRING)
	    call synphoterr (badfile, SW_STR(rg2))

	# Operation performed depends on input type

	if (SW_TYPE(rg1) == SW_LIST && SW_TYPE(rg2) == SW_LIST) {
	    # Combine two lists using union if additive operator
	    # or intersection if multiplicitive operator

	    switch (op) {
	    case Y_ADD, Y_SUB:
		rg3 = syw_union (rg1, rg2)
	    case Y_MUL, Y_DIV:
		rg3 = syw_intersect (rg1, rg2)
	    }

	} else if (SW_TYPE(rg1) == SW_LIST) {
	    # Constant operating on range leaves original range

	    rg3 = syw_copy (rg1)

	} else if (SW_TYPE(rg2) == SW_LIST) {
	    rg3 = syw_copy (rg2)

	} else {
	    # Produce a new constant from math operation on two constants

	    rg3 = syw_copy (rg1)
	    if (IS_INDEFR(SW_VAL(rg1)) || IS_INDEFR(SW_VAL(rg2))) {
		SW_VAL(rg3) = INDEFR

	    } else {
		switch (op) {
		case Y_ADD:
		    SW_VAL(rg3) = SW_VAL(rg1) + SW_VAL(rg2)
		case Y_SUB:
		    SW_VAL(rg3) = SW_VAL(rg1) - SW_VAL(rg2)
		case Y_MUL:
		    SW_VAL(rg3) = SW_VAL(rg1) * SW_VAL(rg2)
		case Y_DIV:
		    SW_VAL(rg3) = SW_VAL(rg1) / SW_VAL(rg2)
		}
	    }
	}

	return (rg3)
end

# SYW_MODE -- Contruct a mode string from an array of descriptors

procedure syw_mode (stack, narg, mode, maxch)

pointer	stack[ARB]	# i: array of descriptors
int	narg		# i: number of descriptors
char	mode[ARB]	# o: mode string
int	maxch		# i: max length of mode string
#--
int	ic, iarg

string	badtype  "Incorrect argument type for function"

int	gstrcpy()

begin
	ic = 1
	do iarg = 1, narg {
	    if (SW_TYPE(stack[iarg]) != SW_STRING)
		call synphoterr (badtype, FUNCNAME)

	    ic = ic + gstrcpy (SW_STR(stack[iarg]), mode[ic], maxch-ic)
	    mode[ic] = ','
	    ic = ic + 1
	}

	if (ic > 1)
	    mode[ic-1] = EOS
	    
end


# SYW_NEGATE -- Implement unary minus operation

pointer procedure syw_negate (rg1)

pointer	rg1		# wavelength set descriptor
#--
pointer	rg2

string	badfile  "Cannot read table"

pointer	syw_copy()

begin
	# A string value indicates that a file was not readable

	if (SW_TYPE(rg1) == SW_STRING)
	    call synphoterr (badfile, SW_STR(rg1))

	# Copy descriptor and negate value if it is a scalar

	rg2 = syw_copy (rg1)
	if (SW_TYPE(rg2) == SW_NUMBER)
	    if (! IS_INDEFR (SW_VAL(rg2)))
		SW_VAL(rg2) = - SW_VAL(rg1)

	return (rg2)
end

# SYW_NEW -- Create a new, empty range set

pointer procedure syw_new (length, type)
int	length		# i: the length of the array
int	type		# i: the type of the array
#--
int	rg

begin
	call calloc (rg, LEN_RGSTRUCT, TY_INT)

	SW_LENGTH(rg) = length
	SW_TYPE(rg) = type

	if (length > 0) {
	    switch (type) {
	    case SW_STRING:
		call malloc (SW_LOPTR(rg), length, TY_CHAR)

	    case SW_NUMBER:
		call malloc (SW_LOPTR(rg), length, TY_REAL)

	    case SW_LIST:
		call malloc (SW_LOPTR(rg), length, TY_REAL)
		call malloc (SW_HIPTR(rg), length, TY_REAL)
		call malloc (SW_KPTR(rg), length, TY_REAL)
	    }
	}

	return (rg)
end


# SYW_NWAVE -- Compute the number of wavelengths in a range interval

int procedure syw_nwave (logspace, lowave, hiwave, minwave, maxwave, nwave)

bool	logspace	# i: wavelengths are spaced logarithmically
real	lowave		# i: start of range interval
real	hiwave		# i: end of range interval
real	minwave		# i: start of total range
real	maxwave		# i: end of total range
int	nwave		# i: number of wavelengths in total range
#--
real	count
int	mwave

begin
	if (logspace) {
	    count = (log10(hiwave)  - log10(lowave)) * 
		    nwave / (log10(maxwave) - log10(minwave))
	} else {
	    count = (hiwave - lowave) * nwave / (maxwave - minwave)
	}

	count = min (count, (hiwave - lowave) / (MIN_DELTA * hiwave))
	mwave = int (count + 0.5)

	return (mwave)
end

# SYW_RANGE -- Retrieve range endpoint given its integer encoding

real procedure syw_range (rg, ix) 

pointer	rg		# i: range descriptor
int	ix		# i: integer encoding
#--
real	range

begin
	if (ix > 0) {
	    range = SW_LORANGE(rg,ix-1)
	} else {
	    range = SW_HIRANGE(rg,-ix-1)
	}

	return	(range)
end

# SYW_REAL -- Convert a string to a floating point number

real procedure syw_real (strbuf)

char	strbuf[ARB]	# i: string to be converted
#--
int	ic
real	rval

string	badnumber "Expected number, got string in synphot expression"

int	ctor()

begin
	ic = 1
	if (ctor (strbuf, ic, rval) == 0)
	    call synphoterr (badnumber, strbuf)

	return rval
end

# SYW_SORT -- Sort a set of ranges from smallest to largest wavelength

procedure syw_sort (rg)

pointer	rg		# u: range set descriptor
#--
int	ir, jr, kr, itemp
real	rtemp

string	badtype  "syw_sort: can only sort lists"

begin
	# Check descriptor type

	if (SW_TYPE(rg) != SW_LIST)
	    call synphoterr (badtype, FUNCNAME)

	# Selection sort

	do ir = 0, SW_LENGTH(rg)-1 {
	    # Find smallest element

	    kr = ir
	    do jr = ir+1, SW_LENGTH(rg)-1 {
		if (SW_LORANGE(rg,jr) < SW_LORANGE(rg,kr) ||
		    (SW_LORANGE(rg,jr) == SW_LORANGE(rg,kr) &&
		     SW_HIRANGE(rg,jr) < SW_HIRANGE(rg,kr)))

		    kr = jr
	    }

	    # Swap smallest element into first location

	    if (kr != ir) {
		rtemp = SW_LORANGE(rg,ir)
		SW_LORANGE(rg,ir) = SW_LORANGE(rg,kr)
		SW_LORANGE(rg,kr) = rtemp

		rtemp = SW_HIRANGE(rg,ir)
		SW_HIRANGE(rg,ir) = SW_HIRANGE(rg,kr)
		SW_HIRANGE(rg,kr) = rtemp

		itemp = SW_KEEPER(rg,ir)
		SW_KEEPER(rg,ir) = SW_KEEPER(rg,kr)
		SW_KEEPER(rg,kr) = itemp
	    }
	}
end


# SYW_STR -- Read a string out of the pseudocode

procedure syw_str (pcode, icode, strbuf, maxch)

int	pcode[ARB]	# i: pseudocode array
int	icode		# u: pointer to current code
char	strbuf[ARB]	# o: output string
int	maxch		# i: maximum length of ouput string
#--
int	ic

begin
	ic = 1
	while (pcode[icode] != EOS) {
	    if (ic <= maxch) {
		strbuf[ic] = pcode[icode]
		ic = ic + 1
	    }

	    icode = icode + 1
	}

	strbuf[ic] = EOS
end

# SYW_TABLE -- Compute the wavelength range of a tabulated passband or spectrum

pointer procedure syw_table (fname)

char	fname[ARB]	# i: table name
#--
int	status
pointer	rg
real	minwave, maxwave

pointer	syw_create(), syw_error()
errchk	tabrange

begin
	call tabrange (fname, minwave, maxwave, status)
	if (status == OK)
	    rg = syw_create (minwave, maxwave)
	else
	    rg = syw_error ()

	return (rg)
end

# SYW_UNION -- Merge two ranges into a new range

pointer procedure syw_union (rg1, rg2)

pointer	rg1	# i: first descriptor
pointer	rg2	# i: second descriptor
#--
int	nr1, nr2, nrange
pointer	rg3

string	badtype  "syw_union: can only take union of lists"

pointer	syw_new()

begin
	# Check descriptors to make sure they are wavelength sets

	if (SW_TYPE(rg1) != SW_LIST || SW_TYPE(rg2) != SW_LIST)
	    call synphoterr (badtype, FUNCNAME)

	# Create new wavelength descriptor

	nrange = SW_LENGTH(rg1) + SW_LENGTH(rg2)
	rg3 = syw_new (nrange, SW_LIST)
	
	# Combine the ranges from the two descriptors

	nr1 = SW_LENGTH(rg1)
	nr2 = SW_LENGTH(rg2)

	call amovr (SW_LORANGE(rg1,0), SW_LORANGE(rg3,0), nr1)
	call amovr (SW_HIRANGE(rg1,0), SW_HIRANGE(rg3,0), nr1)
	call amovi (SW_KEEPER(rg1,0), SW_KEEPER(rg3,0), nr1)

	call amovr (SW_LORANGE(rg2,0), SW_LORANGE(rg3,nr1), nr2)
	call amovr (SW_HIRANGE(rg2,0), SW_HIRANGE(rg3,nr1), nr2)
	call amovi (SW_KEEPER(rg2,0), SW_KEEPER(rg3,nr1), nr2)

	return (rg3)
end

# SYW_WAVE -- Compute a wavelength set for a given set of ranges

procedure syw_wave (rg, wave, nwave)

pointer	rg		# i: range set descriptor
pointer	wave		# o: wavelength set
int	nwave		# u: number of wavelengths in wavelength set
#--
bool	logspace
int	stack[MAXSTACK]
int	ir, jr, nx, ix, jx, top, maxlen, mwave, iwave
pointer	sp, index, npoint
real	minwave, maxwave, rlo, rhi

data	logspace  / true /
string	badexpr  "Synthetic photometry expression too complex"

extern	syw_compare
int	syw_nwave()
real	syw_range()

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (index, 2*SW_LENGTH(rg), TY_INT)
	call salloc (npoint, SW_LENGTH(rg), TY_INT)

	# Sort the wavelength end points

	jr = 0
	do ir = 1, SW_LENGTH(rg) {
	    Memi[index+jr] = ir
	    jr = jr + 1

	    Memi[index+jr] = - ir
	    jr = jr + 1
	}

	nx = jr
	call gqsort (Memi[index], nx, syw_compare, rg)

	# Get minumum and maximum wavelengths

	minwave = syw_range (rg, Memi[index])
	maxwave = syw_range (rg, Memi[index+nx-1])

	# Compute number of wavelengths in each range

	maxlen = nwave
	do ir = 0, SW_LENGTH(rg)-1 {
	    mwave = syw_nwave(logspace, SW_LORANGE(rg,ir), SW_HIRANGE(rg,ir),
			      minwave, maxwave, nwave)

	    if (mwave < MINLENGTH) {
		mwave = syw_nwave (logspace, SW_LORANGE(rg,ir), 
				   SW_HIRANGE(rg,ir), SW_LORANGE(rg,ir), 
				   SW_HIRANGE(rg,ir), MINLENGTH)
		maxlen = maxlen + mwave
	    }

	    Memi[npoint+ir] = mwave
	}

	# Compute non-overlapping wavelength ranges and their wavelength sets

	top = 0
	iwave = 0
	call malloc (wave, maxlen, TY_REAL)

	do ix = 0, nx-1 {
	    rhi = syw_range (rg, Memi[index+ix])
	    ir = abs(Memi[index+ix]) - 1

	    if (ix == 0) {
		rlo = rhi

	    } else if (top == 0) {
		mwave = syw_nwave (logspace, rlo, rhi, minwave, maxwave, nwave)
		call syw_addwave (logspace, rlo, rhi, mwave, 
				  maxlen, iwave, wave)
		rlo = rhi
		
	    } else if (rlo < rhi) {
		jr = stack[top]
		mwave = syw_nwave (logspace, rlo, rhi, SW_LORANGE(rg,jr), 
				   SW_HIRANGE(rg,jr), Memi[npoint+jr])

		call syw_addwave (logspace, rlo, rhi, mwave, 
				  maxlen, iwave, wave)
		rlo = rhi
	    }

	    if (Memi[index+ix] > 0) {
		top = top + 1
		if (top > MAXSTACK)
		    call synphoterr (badexpr, FUNCNAME)

		stack[top] = ir

	    } else {
		for (jx = 1; jx <= top; jx = jx + 1) {
		    if (stack[jx] == ir) 
			break
		}

		for (jx = jx + 1; jx <= top; jx = jx + 1)
		    stack[jx-1] = stack[jx]

		top = top - 1
	    }
	}

	nwave = iwave
	call sfree (sp)
end

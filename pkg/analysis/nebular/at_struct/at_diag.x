# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<evvexpr.h>
include	"../at.h"
include	"../neberr.h"

define	AT_FUNCS	"|J|EXTINCT|"
define	J		1
define	EXTINCT		2

define	LEN_ED_EXPR	2			# expression descriptor
define	ED_N_AT		Memi[$1+0]		# no. atomic descriptors
define	ED_AT		Memi[$1+1]		# ptr to atomic data objects

#define	ED_AT		Memi[ED_AT_PTR($1)+$2-1] # individual at data objects

define	DEBUG		NO

#--------------------------------------------------------------------22 May 97--
.help at_diag.x Apr96 nebular/fivel
.ih
NAME
 at_evdiag - Evaluate the diagnostic expressions for all transitions
  ev_ufunc - Evaluate the named user function
  ev_getop - Called by evvexpr to fetch operand 
  jl_fetch - Return emissivity from named transition
.ih
DESCRIPTION
.ih
ROUTINE DETAILS
.ls double procedure at_evdiag (expr, at)
.ls ARGUMENTS
.ls expr[ARB] (input: char)
Transition expression.
.le
.ls at (input: pointer)
Atomic data object.
.le
.le
.le

.ls procedure ev_getop (ed, opname, o)
.ls ARGUMENTS
.ls ed (input: pointer)
pointer to expression descriptor
.le
.ls opname[ARB] (input: char)
operand name
.le
.ls o (output: pointer)
operand pointer
.le
.le
.le

.ls procedure ev_ufunc (ed, fcname, args, nargs, o)
.ls ARGUMENTS
.ls ed (input: ponter)
pointer to expression object
.le
.ls fcname[ARB] (input: char)
name of user function to be executed
.le
.ls args[ARB] (input: pointer)
pointers to function arguments
.le
.ls nargs (input: int)
no. function arguments
.le
.ls o (output: pointer)
operand pointer
.le
.le
.le
.endhelp
#-------------------------------------------------------------------------------
#  AT_EVDIAG -	Evaluate the diagnostic ratio expression for all transitions.  

real procedure at_evdiag (at, expr)

#  Calling arguments
pointer	at			# I: atomic data object
char	expr[ARB]		# I: constraint expression
real	rval			# O: reduced expression

#  Local variables:
pointer	ed			# expression descriptor
pointer	o			# output operand struct
pointer	sp			# top of stack memory

#  Functions used:
extern	ev_ufunc()		# Evaluate the named user function
extern	ev_getop()		# Evaluate the named operand
pointer	evvexpr()		# evaluate string expression
int	locpr()			# address of procedure

errchk	evvexpr

begin
	# Build the transition descriptor.
	call smark (sp)
	call salloc (ed, LEN_ED_EXPR, TY_STRUCT)
	ED_AT(ed) = at

	# Evaluate the expression; check for divzero, etc.
	rval  = 0.
	o = evvexpr (expr, locpr(ev_getop), ed, locpr(ev_ufunc), ed, EV_RNGCHK)
	rval = O_VALR(o)

	call sfree (sp)
	call mfree (o, TY_STRUCT)
	return (rval)
end


#-------------------------------------------------------------------------------
#  EV_GETOP --	Called by evvexpr to fetch operand 

procedure ev_getop (ed, opname, o)

#  Calling arguments
pointer	ed			# I: expression descriptor
char	opname[ARB]		# I: operand name 
pointer	o			# O: operand pointer

#  Functions used:
#bool	streq()			# are two strings equal?

begin
	# Just copy the string into the operand array buffer. 
	O_TYPE(o) = TY_CHAR
	O_LEN(o)  = 2
	call malloc (O_VALP(o), 2, TY_CHAR)
	call strcpy (opname, O_VALC(o), SZ_FNAME)

end


#-------------------------------------------------------------------------------
#  EV_UFUNC --	Called by evvexpr to evaluate special "constraint" functions 

procedure ev_ufunc (ed, fcname, args, nargs, o)

#  Calling arguments
pointer	ed			# I: pointer to expression descriptor
char	fcname[ARB]		# I: name of user function to be executed
pointer	args[ARB]		# I: pointers to function arguments
int	nargs			# I: no. function arguments
pointer	o			# O: operand pointer

#  Local variables:
pointer	at			# atomic data object
real	emissiv			# transition emissivity
char	errmsg[SZ_LINE]		# text of error message
int	func			# function type
int	lower, upper		# lower, upper levels of transition
pointer	sp, buf			# stack memory pointers
int	v_nargs			# number of expected arguments 

#  Functions used:
double	jl_fetch()		# return transition emissivity
int	strdic()		# return entry in string dictionary
bool	strne()			# are two strings unequal?

begin
	call smark (sp)
	call salloc (buf, SZ_FNAME, TY_CHAR)
	call flush (STDOUT)

	# Look up function in dictionary. 
	call strupr (fcname)
        func = strdic (fcname, Memc[buf], SZ_FNAME, AT_FUNCS)
        if (func > 0 && strne (fcname, Memc[buf]))
            func = 0

	# Abort if the function is not known.
        if (func <= 0) {
	    call sprintf (errmsg, SZ_LINE, "Unknown function `%s'")
		call pargstr (fcname)
	    call error (USER_FUNCTION, errmsg)
	}

	# Verify the correct number of arguments.
	v_nargs = 0
	switch (func) {
	case J:
	    v_nargs = 2
	case EXTINCT:
	    v_nargs = 1
	}
        if (v_nargs > 0 && nargs != v_nargs) {
	    call sprintf (errmsg, SZ_LINE, 
			"Function `J' requires 2 arguments'")
	    call error (USER_FUNCTION, errmsg)
	}

	# Evaluate the functions. 
	switch (func) {
	case J:
	    if (O_TYPE(args[1]) != TY_INT || O_TYPE(args[2]) != TY_INT)
		call error (USER_FUNCTION, 
			"Arguments to J function have incorrect type")

	    # Fetch the emissivity of the referenced transition.
	    upper = O_VALI(args[1])
	    lower = O_VALI(args[2])
	    at = ED_AT(ed,1)

	    iferr (emissiv = jl_fetch (EMISS(at), AT_NLVL(at), upper, lower)) {
	    	call sprintf (errmsg, SZ_LINE, 
		"Expression references invalid transition: %d->%d")
		    call pargi (upper)
		    call pargi (lower)
		call error (USER_FUNCTION, errmsg)
	    }

	    O_VALR(o) = emissiv
	    O_TYPE(o) = TY_REAL
	    O_LEN(o)  = 0			# scalar result

	    if (DEBUG == YES) {
		call printf ("Atomic data descriptor is: %d \n")
		    call pargi (at)
		call printf ("Upper level is: %d \n")
		    call pargi (upper)
		call printf ("Lower level is: %d \n")
		    call pargi (lower)
	    }

	# Other functions not yet supported.
	case EXTINCT:
	    O_TYPE(o) = ERR
	}

	call sfree (sp)
end


#-------------------------------------------------------------------------------
#  JL_FETCH --	Return emissivity from named transition. 

double procedure jl_fetch (jl, sz_jl, upper, lower)

#  Calling arguments
double	jl[sz_jl, sz_jl]	# I: emissivity matrix
int	sz_jl			# I: size of emissivity matrix
int	upper			# I: upper level of transition 
int	lower			# I: lower level of transition
double	emissiv			# O: emissivity

begin
	emissiv = INDEFR

	# Enforce transition selection rules. 
	# Only emission is considered. 
	if (lower >= upper)
	    call error (1, "Transition not allowed")

	# No transitions below ground state, or above highest level modelled.
	else if (lower < 1 || upper > sz_jl)
	    call error (1, "Transition not allowed")

	return (jl[upper, lower])
end



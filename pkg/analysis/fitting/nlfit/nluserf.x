include <error.h>
include	"nlfit.h"

# NL_COMPILE -- Compiles user function.
#
# Uses routine vex_compile in st4gem$lib/xtools/ to parse a fortran 
# expression. The output pointer can be used by nl_userf in latter 
# calls. Upon termination, everything must be cleared by a call to 
# nl_freec.

procedure nl_compile (nl, funct)

pointer	nl		# i: NLFIT structure pointer
char	funct[ARB]	# i: function string

#--
include "nluserf.com"

pointer	vex_compile()

errchk  vex_compile, malloc

begin
	# create common storage
	call malloc (lx, NL_NPTS(nl), TY_REAL)
	call malloc (lcoef, NUPAR, TY_REAL)

	# compile the expression
	NL_UCODE(nl) = vex_compile (funct)

	# enable variable names syntax checking in getval
	optimize = false
end


# NL_FREEC -- Frees areas used by the user function.

procedure nl_freec (nl)

pointer	nl
include "nluserf.com"

errchk  vex_free, mfree

begin
	call mfree (lx, TY_REAL)
	call mfree (lcoef, TY_REAL)
	call vex_free (NL_UCODE(nl))
end


# NL_USERF -- Evaluates pre-compiled user function.
#
# Uses routine vex_eval in st4gem$lib/xtools/ to evaluate a fortran
# expression previously compiled by vex_compile. 

procedure nl_userf (nl, x, z, ydata, ndata, coef, ncoef)

pointer	nl		# i: NLFIT structure pointer.
real	x[ARB]		# i: Indep. variable array.
real	z[ARB]		# o: Computed user function.
real	ydata[ARB]	# i: Dependent variable array.
int	ndata		# i: Number of data points.
real	coef[ARB]	# i: Coefficient array.
int	ncoef		# i: Number of coefficients.

#--
include "nluserf.com"

extern	getvar
int	type

errchk  vex_eval, amovr

begin
	# initialize common before evaluation
	call amovr (x, Memr[lx], ndata)
	call amovr (coef, Memr[lcoef], ncoef)
	lndata = ndata
	lncoef = ncoef

	# evaluate
	iferr (call vex_eval (NL_UCODE(nl), getvar, 0.1d0, type)) {
	    call erract (EA_FATAL)
	    return
	}
	call vex_copyr (NL_UCODE(nl), 0.0, z, ndata)

	# disable variable name syntax checking in following 
	# calls to getvar	
	optimize = true
end


# GETVAR -- Routine called by vex_eval to return the vector associated with
# a variable. Only two possible variable names are accepted, all others
# generate an error abort:
#
#	x   -  The independent variable vector. 
#	ci  -  The i-th function coefficient, where 1 <= i <= ncoef.
#	       In this case, the coefficient value is replicated all over
#	       the buffer elements.
#
# Two versions of the same code are provided, one with variable names
# syntax checking, and another without. The choice is dictated by
# common variable 'optimize'.

procedure getvar (stack, name)

pointer	stack		# u: Stack descriptor
char	name[ARB]	# i: Variable name.
#--
include "nluserf.com"
int	i, j, ival
pointer	buffer

int	strcmp(), ctoi()
pointer	stk_alloc()

string	nameerr	"Error in variable name."

errchk  amovr, amovkr

begin
	if (!optimize) {

	    # syntax checking enabled

	    # independent variable
	    if (strcmp(name, "x") == 0) {
		buffer = stk_alloc (stack, lndata, TY_REAL)
	        call amovr (Memr[lx], Memr[buffer], lndata)

	    # function coefficient
	    } else if ((name[1] == 'c') || (name[1] == 'C')) {
	        i = 2
	        if (ctoi (name, i, ival) == 0)
	            call error (0, nameerr)
	        if ((ival <= lncoef) && (ival > 0)) {
	            if (Memr[lcoef+ival-1] != INDEF) {
			buffer = stk_alloc (stack, lndata, TY_REAL)
	                call amovkr (Memr[lcoef+ival-1], Memr[buffer], lndata)
	            } else {
	                call error (0, nameerr)
	            }
	        } else 
	            call error (0, nameerr)
	    } else
	        call error (0, nameerr)
	} else {

	    # syntax checking disabled

	    # independent variable
	    if (name[1] == 'x') {
		buffer = stk_alloc (stack, lndata, TY_REAL)
	        call amovr (Memr[lx], Memr[buffer], lndata)

	    # function coefficient
	    } else {
	        i = 2
	        j = ctoi (name, i, ival)
	        ival = ival - 1
		buffer = stk_alloc (stack, lndata, TY_REAL)
	        call amovkr (Memr[lcoef+ival], Memr[buffer], lndata)
	    }
	}
end

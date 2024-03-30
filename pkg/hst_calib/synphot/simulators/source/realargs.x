include "function.h"

#* HISTORY *
#* B.Simon	17-Feb-95	original
#* B.Simon	04-May-95	Changed function descriptor

# REALARGS -- Copy and verify function arguments to a real array

int procedure realargs (func, offset, arglist, nlist)

pointer	func		# i: object function descriptor
int	offset		# i: first parameter to copy
real	arglist[ARB]	# o: array of real arguments
int	nlist		# i: number of arguments to copy
#--
int	ilist, jlist
int	btoi()

begin
	jlist = offset
	for (ilist = 1; ilist <= nlist; ilist = ilist + 1) {
	    if (jlist > FUN_NPAR(func))
		break

	    if (FUN_TYPE(func,jlist) != TY_REAL)
		break

	    arglist[ilist] = FUN_NUM(func,jlist)
	    jlist = jlist + 1
	}

	return (btoi (ilist > nlist && jlist > FUN_NPAR(func)))
end

include "function.h"

#* HISTORY *
#* B.Simon	04-May-95	Original

# NAMEFUNC -- Recreate the name string from the function descriptor

procedure namefunc (func, name, maxch)

pointer	func		# i: function descriptor
char	name[ARB]	# o: function name string
int	maxch		# i: max length name string
#--
int	ic, ipar
pointer	sp, value

int	gstrcpy()

begin
	# Allocate memory for value string

	call smark (sp)
	call salloc (value, SZ_FNAME, TY_CHAR)

	# Write function name

	ic =  1 + gstrcpy (FUN_STR(func,1), name, maxch)

	# Write parenthesized list of function arguments

	ic = ic + gstrcpy ("(", name[ic], maxch-ic)

	do ipar = 2, FUN_NPAR(func) {
	    if (ipar > 2)
		ic = ic + gstrcpy (",", name[ic], maxch-ic)

	    if (FUN_TYPE(func,ipar) == TY_CHAR) {
		ic = ic + gstrcpy (FUN_STR(func,ipar), name[ic], maxch-ic)

	    } else {
		call sprintf (Memc[value], SZ_FNAME, "%g")
		call pargr (FUN_NUM(func, ipar))

		ic = ic + gstrcpy (Memc[value], name[ic], maxch-ic)
	    }
	}

	ic = ic + gstrcpy (")", name[ic], maxch-ic)

	call sfree (sp)
end

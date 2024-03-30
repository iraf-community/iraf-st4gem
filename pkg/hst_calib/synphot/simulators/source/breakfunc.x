include	"function.h"

#* HISTORY *
#* B.Simon	17-Feb-95	original
#* B.Simon	04-May-95	Changed function descriptor

# BREAKFUNC -- Parse a function into its name and parameters

pointer procedure breakfunc (str)

char	str[ARB]	# i: function string
#--
char	lpar, rpar
int	ic, jc, nc, lpos, rpos, nparam, ipar
pointer	sp, arg, arglist, func, param, partype
real	value

data	lpar, rpar  / '(', ')' /

int	stridx(), strlen(), ctor(), word_count(), word_fetch()

begin
	# Allocate memory for strings

	call smark (sp)
	call salloc (arg, SZ_FNAME, TY_CHAR)
	call salloc (arglist, SZ_LINE, TY_CHAR)

	# Separate function name and argument list

	lpos = stridx (lpar, str)
	rpos = stridx (rpar, str)

	if (rpos <= lpos)
	    rpos = strlen (str) + 1

	nc = (rpos - lpos) - 1
	call strcpy (str[lpos+1], Memc[arglist], nc)

	# Allocate memory for function structure

	nparam = 1 + word_count (Memc[arglist])

	call malloc (func, LEN_FUNSTRUCT, TY_STRUCT)
	call malloc (param, nparam, TY_INT)
	call malloc (partype, nparam, TY_INT)

	FUN_NPAR(func) = nparam
	FUN_TYPARY(func) = partype
	FUN_PARARY(func) = param

	# The first parameter is the function name

	Memi[partype] = TY_CHAR
	call malloc (FUN_PARAM(func,1), SZ_FNAME, TY_CHAR)

	if (lpos == 0) {
	    FUN_STR(func,1) = EOS
	} else {
	    call strcpy (str, FUN_STR(func,1), lpos-1)
	}

	# Extract remaining parameters from argument list

	ic = 1
	ipar = 2
	while (word_fetch (Memc[arglist], ic, Memc[arg], SZ_LINE) > 0) {

	    jc = 1
	    nc = ctor (Memc[arg], jc, value)

	    if (Memc[arg+jc-1] != EOS) {
		FUN_TYPE(func,ipar) = TY_CHAR
		call malloc (FUN_PARAM(func,ipar), SZ_FNAME, TY_CHAR)
		call strcpy (Memc[arg], FUN_STR(func,ipar), SZ_FNAME)

	    } else {
		FUN_TYPE(func,ipar) = TY_REAL
		call malloc (FUN_PARAM(func,ipar), 1, TY_REAL)
		FUN_NUM(func,ipar) = value
	    }

	    ipar = ipar + 1
	}

	return (func)
end

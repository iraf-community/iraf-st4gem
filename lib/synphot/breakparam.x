include "libsynphot.h"

# BREAKPARAM -- Extract parameters from parameterized string

procedure breakparam (keyword, mxparam, nparam, param)

char	keyword[ARB]	# u: Instrument mode keyword
int	mxparam		# i: Maximum number of parameters in string
int	nparam		# o: Number of parameters in string
real	param[ARB]	# o: Parameters extracted from keyword
#--
char	pchar
int	ic, jc, nc

data	pchar	/ PCH /

string	overparam  "Too many parameters in string"
string  badparam   "Illegal character in parameter"
string  noparam    "Parameter value is missing"

int	stridx(), ctor()
errchk	synphoterr

begin
	ic = stridx (pchar, keyword)
	if (ic == 0) {
	    nparam = 0

	} else {
	    jc = ic
	    nparam = 0

	    repeat {
		if (nparam == mxparam)
		    call synphoterr (overparam, keyword)

		# Add a trailing char to keyword for each parameter

		keyword[ic] = pchar
		nparam = nparam + 1
		ic = ic + 1
		jc = jc + 1

		# Convert parameter to real number and add to list

		nc = ctor (keyword, jc, param[nparam])
		if (nc == 0)
		    call synphoterr (noparam, keyword)

	    } until (keyword[jc] != pchar)

	    keyword[ic] = EOS

	    # Error exit if non-numeric character in parameter

	    if (keyword[jc] != EOS)
		call synphoterr (badparam, keyword)
	}

end

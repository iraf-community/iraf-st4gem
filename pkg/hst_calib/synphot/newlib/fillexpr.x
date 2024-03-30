include <ctype.h>

#* HISTORY *
#* B.Simon	29-Apr-94	original

# FILLEXPR -- Replace the variables in a synphot expression with their values

procedure fillexpr (istr, ostr, maxch)

char	istr[ARB]	# i: expression with variables
char	ostr[ARB]	# o: expression with variables replaced 
int	maxch		# i: maximum length of ostr
#--
bool	inval
int	ic, jc, oc
pointer	sp, value
real	rval

begin
	call smark (sp)
	call salloc (value, SZ_FNAME, TY_CHAR)

	ic = 1
	oc = 1
	inval = false

	# This procedure toggles between two states, depending on inval.
	# When inval is true, characaters are being copied from
	# the variable value to the output string. This state terminates
	# when the EOS is seen at the end of the value string.
	# When inval is false, characters are copied from the input string
	# to the output string. This state terminates when a dollar sign 
	# followed by a digit is seen.

	while (oc <= maxch) {
	    if (inval) {
		if (Memc[value+jc] == EOS) {
		    inval = false

		} else {
		    ostr[oc] = Memc[value+jc]
		    oc = oc + 1
		    jc = jc + 1
		}

	    } else {
		if (istr[ic] == EOS)
		    break

		if (istr[ic] == '$') {
		    if (IS_DIGIT(istr[ic+1])) {
			call getsynvar (TO_INTEG(istr[ic+1]), rval)
			call sprintf (Memc[value], SZ_FNAME, "%g")
			call pargr (rval)

			inval = true
			ic = ic + 2
			jc = 0

		    } else {
			ostr[oc] = istr[ic]
			oc = oc + 1
			ic = ic + 1
		    }

		} else {
		    ostr[oc] = istr[ic]
		    oc = oc + 1
		    ic = ic + 1
		}
	    }
	}

	ostr[oc] = EOS
	call sfree (sp)
end

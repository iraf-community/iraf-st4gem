include	<ctype.h>

# STRJUST -- Justify a text string
#
# Justifying a string means converting to lower case and removing 
# whitespace. The output string may be the same as the input string.
#
# B.Simon	15-Aug-88	First Code

procedure strjust (in, out, maxch)

char	in[ARB]		# i: Input string
char	out[ARB]	# o: Output string
int	maxch		# i: Length of output string
#--
char	ch
int	ic, oc

begin
	oc = 1
	for (ic = 1; in[ic] != EOS; ic = ic + 1) {
	    ch = in[ic]
	    if (! IS_WHITE(ch)) {
		if (IS_UPPER(ch))
		    out[oc] = TO_LOWER(ch)
		else
		    out[oc] = ch
		oc = oc + 1
		if (oc > maxch)
		    break
	    }
	}

	out[oc] = EOS
end

include	<ctype.h>

# STRFIX -- Remove whitspace from a string and convert to lower case

procedure strfix (str)

char	str[ARB]	# u: string to convert
#--
int ic, jc

begin
	jc = 1
	for (ic = 1; str[ic] != EOS; ic = ic + 1) {
	    if (IS_WHITE(str[ic]))
		next

	    if (IS_UPPER(str[ic])) {
		str[jc] = TO_LOWER(str[ic])

	    } else if (jc < ic) {
		str[jc] = str[ic]
	    }

	    jc = jc + 1
	}

	str[jc] = EOS
end


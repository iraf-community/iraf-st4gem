include <ctype.h>

#* HISTORY *
#* B.Simon	13-Apr-95	Copied from datefield in cdbs

# TIMEFIELD -- Get the next field from a time string

procedure timefield (time, ic, value, field, maxch)

char	time[ARB]	# i: time string
int	ic		# u: current charcater in field
int	value		# o: integer value of field, if any
char	field[ARB]	# o: output field
int	maxch		# i: declared length of output field
#--
int	jc

begin
	value = 0
	field[1] = EOS

	while (! IS_ALNUM(time[ic]) && time[ic] != EOS)
	    ic = ic + 1


	for (jc = 1; jc <= maxch && IS_ALNUM(time[ic]); jc = jc + 1) {
	    field[jc] = time[ic]

	    if (IS_DIGIT(time[ic])) {
		if (!IS_INDEFI (value))
		    value = 10 * value + TO_INTEG(time[ic])
	    } else {
		value = INDEFI
	    }

	    ic = ic + 1
	}

	field[jc] = EOS
end

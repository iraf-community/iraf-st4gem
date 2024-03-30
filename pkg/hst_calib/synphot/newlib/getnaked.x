#* HISTORY *
#* B.Simon	13-Oct-94	original

# GETNAKED -- Extract the naked obsmode expression from a synphot expression

procedure getnaked (input, output, maxch)

char	input[ARB]	# i: input obsmode expression
char	output[ARB]	# o: output obsmode expression
int	maxch		# i: maximum length of output
#--
int	ic, jc, state, strict
pointer	sp, token

bool	streq()
int	syntok(), gstrcpy()

begin
	# Allocate memory to hold token

	call smark (sp)
	call salloc (token, maxch, TY_CHAR)

	# Fetch each token from the obsmode string

	ic = 1
	jc = 1
	state = 0
	strict = NO
	while (syntok (input, ic, strict, Memc[token], maxch) > 0) {
	    call strlwr (Memc[token])

	    # State machine to skip band() tokens

	    switch (state) {
	    case 0:
		if (streq (Memc[token], "band")) {
		    state = 2
		    next
		}
	    case 1:
		if (streq (Memc[token], ")")) {
		    state = 0
		    next
		}
	    case 2:
		if (streq (Memc[token], "(")) {
		    state = 1
		    next
		}
	    }

	    # Copy other tokens  to output. Add a space 
	    # for the sake of space separated tokens

	    jc = jc + gstrcpy (Memc[token], output[jc], maxch-jc)
	    output[jc] = ' '
	    jc = jc + 1
	}

	output[jc] = EOS
	call sfree (sp)
end

#* HISTORY *
#* B.Simon	25-May-94	original

# EXP_REWRITE -- Rewrite a synphot expression in standard syntax

procedure exp_rewrite (input, output, maxch)

char	input[ARB]	# i: Input expression
char	output[ARB]	# o: Output expression
int	maxch		# i: Maximum length of output expression
#--
bool	blank
int	ic, jc
pointer	sp, word

bool	is_simple()
int	word_fetch(), gstrcpy()

begin
	# Allocate memory for temporary string

	call smark (sp)
	call salloc (word, SZ_FNAME, TY_CHAR)

	# A simple expression is a list of obsmode keywords separated by
	# blanks or commas. The keywords are inserted in the band() 
	# function and the keywords are separated by commas. This is the
	# syntax synphot requries for a passband.
	
	if (is_simple (input)) {
	    jc = gstrcpy ("band(", output, maxch) + 1

	    ic = 1
	    blank = true
	    while (word_fetch (input, ic, Memc[word], SZ_FNAME) > 0) {
		if (jc >= maxch)
		    break

		blank = false

		jc = jc + gstrcpy (Memc[word], output[jc], maxch-jc)
		output[jc] = ','
		jc = jc + 1
	    }

	    if (! blank)
		jc = jc - 1

	    output[jc] = ')'
	    output[jc+1] = EOS

	} else {
	    call strcpy (input, output, maxch)
	}

	call sfree (sp)
end

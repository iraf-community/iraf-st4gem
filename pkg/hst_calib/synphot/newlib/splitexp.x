#* HISTORY *
#* B.Simon	13-May-94	original
#* B.Simon	25-May-94	added call to exp_rewrite

# SPLITEXP -- Split a synphot expression in two at the first minus token

procedure splitexp (command, nexp, exp1, exp2, maxch)

char	command[ARB]	# i: Synphot expression
int	nexp		# o: Number of expressions (1 or 2)
char	exp1[ARB]	# o: First subexpression
char	exp2[ARB]	# o: Second subexpression
int	maxch		# i: Maximum length of subexpressions
#--
int	ic, jc, found, strict
pointer	sp, exp, token

bool	streq()
int	syntok()

begin
	# Allocate memory for temporary string

	call smark (sp)
	call salloc (exp, maxch, TY_CHAR)
	call salloc (token, maxch, TY_CHAR)

	# Replace "naked" obsmodes with calls to the band() function

	call exp_rewrite (command, Memc[exp], maxch)

	# Break command into tokens, look for minus sign
	# The tokenizer has two modes, strict and loose, and
	# sets the value of strict accordingly

	ic = 1
	jc = 1
	found = NO
	strict = YES

	while (syntok (Memc[exp], ic, strict, Memc[token], maxch) > 0) {
	    if (streq (Memc[token], "-")) {
		while (Memc[exp+ic-1] != EOS) {
		    if (Memc[exp+ic-1] > ' ')
			break

		    ic = ic + 1
		}
		    
		found = YES
		break
	    }

	    jc = ic
	}

	# Break input string in two at the minus sign
	# If not found, copy entire string to first output

	if (found == NO) {
	    nexp = 1
	    call strcpy (Memc[exp], exp1, maxch)
	    exp2[1] = EOS

	} else {
	    nexp = 2
	    call strcpy (Memc[exp], exp1, jc-1)
	    call strcpy (Memc[exp+ic-1], exp2, maxch)
	}

	call sfree (sp)

end

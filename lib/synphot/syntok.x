include <ctype.h>
define	BLANK		' '

#* HISTORY *
#* B.Simon	13-May-94	original

# SYNTOK -- Extract the next token from a synphot command

int procedure syntok (command, ic, strict, token, maxch)

char	command[ARB]	# i: Command string
int	ic		# u: Position in command string
int	strict		# u: Strict interpretation of tokens
char	token[ARB]	# o: Token string
int	maxch		# i: Maximum length of token string
#--
bool	digit
int	jc, nc
real	rval

string	loosestr  "(),"
string	strictstr "(),*+-"
string	loosefunc "band,cat,icat,spec,thru"

bool	streq()
int	ctor(), stridx(), strlen(), word_match()

begin
	# Skip over leading white space

	while (command[ic] <= BLANK) {
	    if (command[ic] == EOS)
		break

	    ic = ic + 1
	}

	# Check to see if token is a number

	if (strict == NO) {
	    digit = false

	} else {
	    if (command[ic] == '-') {
		jc = ic + 1
	    } else {
		jc = ic
	    }

	    digit = command[jc] == '.' || IS_DIGIT(command[jc])
	}

	# Use ctor to get the token if it is a number
	# Otherwise, read characters until one is found in the stop set

	if (digit) {
	    jc = 1
	    nc = ctor (command[ic], jc, rval)

	    nc = min (nc, maxch)
	    ic = ic + nc

	    call strcpy (command[ic], token, nc)

	} else {
	    jc = 1
	    while (command[ic] > BLANK) {
		if (strict == NO) {
		    if (stridx (command[ic], loosestr) > 0)
			break
		} else {
		    if (stridx (command[ic], strictstr) > 0)
			break
		}
		if (jc > maxch)
		    break

		token[jc] = command[ic]
		ic = ic + 1
		jc = jc + 1
	    }

	    # A  character in the stop set is a single character token

	    if (jc == 1) {
		token[1] = command[ic]
		ic = ic + 1
		jc = 2
	    }

	    token[jc] = EOS
	    nc = strlen (token)
	}

	# Change the value of strict based on the new token

	if (word_match (token, loosefunc) > 0)
	    strict = NO

	if (streq (token, ")")) 
	    strict = YES

	return (nc)
end


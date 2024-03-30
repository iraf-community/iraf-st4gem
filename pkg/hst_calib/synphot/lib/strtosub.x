# STRTOSUB -- Parse a substring from main string.  Start at ic and parse
# to next occurence of 'ch'.  ic is left pointing to the character after
# the delimiter character ch.  ch is not included in the parsed substring.

int procedure strtosub( str, ip, ch, outstr, maxch)

char	str[ARB]	# i: Input string
int	ip		# i: Pointer to beginning of next substring
char	ch		# i: Character that delimits a word
char	outstr[ARB]	# o: Output substring
int	maxch		# i: Maximum length of str

int	ic, jc, max
int	strlen()

begin

	# Copy input to output while input not equal to delimiter char

	jc = 1
	max = min ( maxch, strlen(str) + 1 )
	for ( ic = ip; ic < max && str[ic] != ch && str[ic] != EOS;
	      ic = ic + 1 ) {
	   outstr[jc] = str[ic]
	   jc = jc + 1
	}
	outstr[jc] = EOS

	ip = ic + 1
	return ( jc - 1 )

end

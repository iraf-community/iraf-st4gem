# SCOUNT -- count the number of occurences of a substring in a string

int procedure scount( str, patstr, maxch )

char	str[ARB]	# i: input string
char	patstr[ARB]	# i: pattern to search for
int	maxch		# i: length of str

int	ic, npat, idx, shift
int	strlen(), strsearch()
begin

	shift = strlen(patstr)
	ic = 1
	for( npat=0; ic < maxch; npat = npat + 1 ) {
	   idx = strsearch( str[ic], patstr )
	   if ( idx > 0 )
	      ic = ic + idx + shift - 1
	   else
	      return ( npat )
	}

	return (npat)
end

# UNBANG -- Remove a bang at the begining of a string.  Returns unbanged str
#  as outstring and number of chars as procedure value

int procedure unbang(instring, outstring, maxch )

char	instring[ARB]	# input string
char	outstring[ARB]	# output string
int	maxch		# max length of outstring

int	ic, jc
int	stridx()

begin

	jc = 0
	for (ic=stridx("!",instring) + 1; instring[ic] != EOS; ic = ic + 1 )
	   if ( ic <= maxch ) {
	      jc = jc + 1
	      outstring[jc] = instring[ic]
	   }

	outstring[jc+1] = EOS
	return ( jc )
end

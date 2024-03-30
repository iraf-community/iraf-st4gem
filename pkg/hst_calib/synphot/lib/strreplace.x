# STRREPLACE -- Replace one substring of a string with another.  Returns
# number of search strings replaced.  Returns -1 if maxch exceeded.

int procedure strreplace( sstr, rstr, mainstr, maxch )

char	sstr[ARB]	# i: Search substring to be replaced
char	rstr[ARB]	# i: Replacement substring
char	mainstr[ARB]	# io: String to be acted upon
int	maxch		# i: Max length of mainstr

int	ic, idx, nstr, rlen
int	strsearch(), strlen()
char	temp[SZ_LINE]

begin
	
	rlen = strlen( rstr )
	ic = 1
	nstr = 0
	while ( ic < maxch ) {

	   # Find next occurence of search string
	   idx = strsearch( mainstr[ic], sstr )
	   if ( idx > 0 ) {	# Found one

	      nstr = nstr + 1

	      # Point to position after sstr
	      ic = ic + idx - 1

	      # Copy remainder of mainstr to temp for safekeeping and
	      # and terminate mainstr
	      call strcpy( mainstr[ic], temp, maxch)

	      # Terminate mainstr at initial position of sstr
	      mainstr[ic - strlen(sstr)] = EOS

	      # Concatenate rstr
	      call strcat(rstr, mainstr, maxch)

	      # Put back remainder of mainstr
	      call strcat(temp, mainstr, maxch)

	   # sstr not found so return

	   } else
	      return ( nstr )
	}
	return ( -1 )
	   
end

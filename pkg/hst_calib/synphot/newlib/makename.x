#* HISTORY *
#* B.Simon	28-Apr-94	original

# MAKENAME -- Create a unique name from a root and an integer

procedure makename (root, iexp, nexp, name, maxch)

char	root[ARB]	# i: root of name
int	iexp		# i: integer to append
int	nexp		# i: total number of names
char	name[ARB]	# o: output name
int	maxch		# i: maximum length of name
#--

begin
	if (nexp <= 1) {
	    call strcpy (root, name, maxch)
	} else {
	    call sprintf (name, maxch, "%s%d")
	    call pargstr (root)
	    call pargi (iexp)
	}
end


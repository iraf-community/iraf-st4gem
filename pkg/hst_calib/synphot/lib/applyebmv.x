# APPLYEBMV -- Apply extinction to a command string

procedure applyebmv (instring, outstring, ebmv, maxch)

char	instring[ARB]	# i: input command without extinction
char	outstring[ARB]	# o: output command with extinction
real	ebmv		# i: value of extinction
int	maxch		# i: maximum number of characters in outstring
#--
pointer	sp, temp

string	simplefmt  "band(%s)*ebmv(%f)"
string	complexfmt "(%s)*ebmv(%f)"

bool	is_simple()

begin
	# If ebmv is zero, just copy string and return

	if (ebmv == 0.0) {
	    call strcpy( instring, outstring, maxch )
	    return
	}

	# If ebmv is not zero, insert the extinction function

	call smark (sp)
	call salloc (temp, maxch, TY_CHAR)

	call strcpy (instring, Memc[temp], maxch)

	if (is_simple (Memc[temp])) {
	    call sprintf( outstring, maxch, simplefmt )
	    call pargstr( Memc[temp] )
	    call pargr( ebmv )
	} else {
	    call sprintf( outstring, maxch, complexfmt )
	    call pargstr( Memc[temp] )
	    call pargr( ebmv )
	}

	call sfree (sp)
end

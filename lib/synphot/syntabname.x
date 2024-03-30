# SYNTABNAME -- Get the name of a synphot table from it descriptor

procedure syntabname (tp, tabname, maxch)

pointer	tp		# i: table descriptor
char	tabname[ARB]	# o: table name
int	maxch		# i: maximum length of table name
#--
int	ic, nc

int	strlen()

begin
	# Get table name

	call tbtnam (tp, tabname, maxch)

	# Strip off bracketed section

	nc = strlen (tabname)
	if (tabname[nc] == ']') {
	    do ic = nc-1, 1, -1 {
		if (tabname[ic] == '[') {
		    tabname[ic] = EOS
		    break
		}
	    }
	}
end


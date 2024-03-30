define	LBRAK		'['
define	RBRAK		']'
define	BSLASH		'\\'

# BREAKCOMP -- Break filename into table and column names

procedure breakcomp (filename, tabname, colname, maxch)

char	filename[ARB]	# i: file name followed by bracketed column name
char	tabname[ARB]	# o: table name
char	colname[ARB]	# o: column name
int	maxch		# i: max length of output strings
#--
int	ic, leftpos, rightpos, tablen, collen

string	badbraks "Right bracket missing from filename"
string	toolong  "File name too long"

errchk	synphoterr

begin
	# Find positions of brackets which contain column name

	leftpos = 0
	rightpos = 0

	for (ic = 1; filename[ic] != EOS; ic = ic + 1) {
	    if (filename[ic] == BSLASH  && filename[ic+1] != EOS) {
		ic = ic + 1

	    } else if (filename[ic] == LBRAK) {
		leftpos = ic

	    } else if (filename[ic] == RBRAK && leftpos > 0) {
		rightpos = ic
		break
	    }
	}

	# If no brackets, set column name to null string
	# Otherwise, break file name into table and column names

	if (leftpos == 0) {
	    call strcpy (filename, tabname, maxch)
	    colname[1] = EOS

	} else if (rightpos == 0) {
	    call synphoterr (badbraks, filename)

	} else {
	    tablen = leftpos - 1
	    collen = rightpos - (leftpos + 1)

	    if (tablen > maxch || collen > maxch)
		call synphoterr (toolong, filename)

	    call strcpy (filename, tabname, tablen)
	    call strcpy (filename[leftpos+1], colname, collen)
	}

end

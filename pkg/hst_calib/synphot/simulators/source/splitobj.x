define	done_		90

#* HISTORY *
#* B.Simon	17-Feb-95	original

# SPLITOBJ -- Split an input line into its elements

int procedure splitobj (line, x, y, mag, spec, object, maxch)

char	line[ARB]	# u: Line from input file
real	x		# o: X position
real	y		# o: Y position
real	mag		# o: Magnitude
char	spec[ARB]	# o: Object spectrum
char	object[ARB]	# o: Object shape
int	maxch		# i: Length of output strings
#--
char	quote
int	ic, nc, nword

string	badline  "Incomplete input"
string	dots     " ..."

int	ctor(), ctowrd(), strlen()

begin
	# Remove comment from line

	quote = EOS
	for (ic = 1; line[ic] != EOS; ic = ic + 1) {
	    if (quote == EOS) {
		if (line[ic] == '\'' || line[ic] == '\"') {
		    quote = line[ic]
		} else if (line[ic] == '#') {
		    line[ic] = EOS
		    break
		}

	    } else {
		if (quote == line[ic])
		    quote = EOS
	    }
	}

	# Initialize output

	x = INDEFR
	y = INDEFR
	mag = INDEFR
	spec[1] = EOS
	object[1] = EOS

	# Read four fields from line: two numbers and two strings

	ic = 1

	if (ctor (line, ic, x) == 0) {
	    nword = 0
	    goto done_
	}

	if (ctor (line, ic, y) == 0) {
	    nword = 1
	    goto done_
	}

	if (ctor (line, ic, mag) == 0) {
	    nword = 2
	    goto done_
	}

	if (ctowrd (line, ic, spec, maxch) == 0) {
	    nword = 3
	    goto done_
	}

	if (ctowrd (line, ic, object, maxch) == 0) {
	    nword = 4
	    goto done_
	}

	nword = 5

	# Check number of fields read

done_	if (nword > 0 && nword < 3) {
	    nc = strlen (line)
	    if (nc > 64) {
		line[60] = EOS
		call strcat (dots, line, nc)
	    }

	    call printerr_str (badline, line)
	}

	return (nword)
end
			  

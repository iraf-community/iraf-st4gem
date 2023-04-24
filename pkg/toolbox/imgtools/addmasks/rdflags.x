# RD_FLAGS -- Read the list of flags from the input parameter
#
# B.Simon	04-Apr-91	First Code

procedure rd_flags (flags, flist, nlist, maxlist)

char	flags[ARB]	# i: flag parameter
int	flist[ARB]	# o: list of flag values
int	nlist		# o: number of flags in list
int	maxlist		# i: maximum length of list
#--
char	cmt
int	ic, jc, fd
pointer	sp, line, word, errmsg

data	cmt	/ '#' /

string	toobig  "List of flag values exceeds maximum length (%d)"
string	notpos  "Illegal flag value (%s)"

int	word_fetch(), getline(), open(), ctoi(), stridx()

begin
	# Allocate dynamic memory for strings
	call smark (sp)
	call salloc (line, SZ_LINE, TY_CHAR)
	call salloc (word, SZ_LINE, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	# Skip over leading blanks in flags parameter
	for (ic = 1; flags[ic] != EOS; ic = ic + 1) {
	    if (flags[ic] > ' ')
		break
	}

	# Use the line array to disguise differences btw reading
	# flag values directly and from a file
	if (flags[ic] == '@') {
	    fd = open (flags[ic+1], READ_ONLY, TEXT_FILE)
	    Memc[line] = EOS
	} else {
	    fd = 0
	    call strcpy (flags[ic], Memc[line], SZ_LINE)
	}

	# Loop over each word in the flag list
	# Words are separated by blanks or commas
	# Anything following a '#' is ignored
	ic = 1
	nlist = 0

	repeat {
	    # If there are no more words on the line, read next line
	    # or break if there was no file or no more lines in the file
	    if (Memc[line+ic-1] == EOS) {
		if (fd == 0) {
		    break
		} else if (getline (fd, Memc[line]) == EOF) {
		    break
		} else {
		    ic = 1
		    # Chop line at the comment character
		    jc = stridx (cmt, Memc[line])
		    if (jc > 0)
			Memc[line+jc-1] = EOS
		}
	    }

	    # Get next word from the line
	    if (word_fetch (Memc[line], ic, Memc[word], SZ_LINE) > 0) {
		# Increment number of flags, check for overflow
		nlist = nlist + 1
		if (nlist > maxlist) {
		    call sprintf (Memc[errmsg], SZ_LINE, toobig)
		    call pargi (maxlist)

		    call error (1, Memc[errmsg])
		}

		# Convert to integer, check for illegal values
		jc = 1
		jc = ctoi (Memc[word], jc, flist[nlist])
		if (flist[nlist] <= 0) {
		    call sprintf (Memc[errmsg], SZ_LINE, notpos)
		    call pargstr (Memc[word])

		    call error (1, Memc[errmsg])
		}
	    }
	}

	if (fd != 0)
	    call close (fd)

	call sfree (sp)
end

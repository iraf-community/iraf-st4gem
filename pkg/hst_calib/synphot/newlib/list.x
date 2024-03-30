define	LEN_LSTRUCT	2

define	L_FILE		Memi[$1]	# Spool file holding list
define	L_SIZE		Memi[$1+1]	# Number of lines in spool file


#* HISTORY *
#* B.Simon	28-Apr-94	original
#* B.Simon	11-May-94	only allow whole line comments

# FREELIST -- Close the list and free associated memory

procedure freelist (list)

pointer	list		# i: list descriptor
#--

begin
	call close (L_FILE(list))
	call mfree (list, TY_STRUCT)
end


# NUMLIST -- Get the number of lines in the list

int procedure numlist (list)

pointer	list		# i: list descriptor
#--

begin
	return (L_SIZE(list))
end

# NXTLIST -- Get the next line from the list

int procedure nxtlist (list, str, maxch)

pointer	list		# i: list descriptor
char	str[ARB]	# o: line from list
int	maxch		# i: maximum length of string
#--
bool	blank
int	status
pointer	sp, line,ch

int	getline()

begin
	# Allocate string to hold line

	call smark (sp)
	call salloc (line, SZ_LINE, TY_CHAR)

	# Read file until non-blank, non-comment line found

	blank = true
	while (blank) {
	    if (getline (L_FILE(list), Memc[line]) == EOF)
		break

	    # Remove trailing newline and comments

	    for (ch = line; Memc[ch] != EOS; ch = ch + 1) {
		if (Memc[ch] == '\n') {	
		    Memc[ch] = EOS
		    break
		} else if (blank && Memc[ch] == '#') {
		    Memc[ch] = EOS
		    break
		} else if (Memc[ch] > ' ') {
		    blank = false
		}
	    }
	}

	# Copy line to output

	if (blank) {
	    call seek (L_FILE(list), BOF)
	    str[1] = EOS
	    status = EOF
	} else {
	    call strcpy (Memc[line], str, maxch)
	    status = OK
	}

	call sfree (sp)
	return (status)
end

# RDLIST -- Read a one column list from an ascii file

pointer procedure rdlist (str)

char	str[ARB]	# i: string or filename
#--
char	nl
int	ic, fd, ispool
pointer	sp, spool, line, list

data	nl	/ '\n' /
data	ispool	/ 1 /

int	open(), getline()

begin
	call smark (sp)
	call salloc (spool, SZ_FNAME, TY_CHAR)
	call salloc (line, SZ_LINE, TY_CHAR)

	# Allocate list structure

	call malloc (list, LEN_LSTRUCT, TY_STRUCT)

	call sprintf (Memc[spool], SZ_FNAME, "spool_%d")
	call pargi (ispool)
	ispool = ispool + 1

	L_FILE(list) = open (Memc[spool], READ_WRITE, SPOOL_FILE)
	L_SIZE(list) = 0

	# Search for first non-white character in string

	for (ic = 1; str[ic] != EOS; ic = ic + 1)
	    if (str[ic] > ' ')
		break

	# Read list file into spool file

	if (str[ic] != '@') {
	    # Handle case where string is not filename

	    call putline (L_FILE(list), str[ic])
	    call putc (L_FILE(list), nl)
	    L_SIZE(list) = L_SIZE(list) + 1

	} else {
	    # Read each line of file individually

	    fd = open (str[ic+1], READ_ONLY, TEXT_FILE)
	    while (getline (fd, Memc[line]) != EOF) {
		if (Memc[line] != '#' && Memc[line] != '\n') {
		    call putline (L_FILE(list), Memc[line])
		    L_SIZE(list) = L_SIZE(list) + 1
		}
	    }

	    call close (fd)
	}

	# Rewind spool file

	call flush (L_FILE(list))
	call seek (L_FILE(list), BOF)

	call sfree (sp)
	return (list)
end

#* HISTORY *
#* B. Simon	05-Aug-94	original

# WRTCOMMAND -- Write the command string piecewise to a file

procedure wrtcommand (putfunc, fd, maxch, command)

extern	putfunc		# i: function to write string
int	fd		# i: file descriptor
int	maxch		# i: maximum length of string piece
char	command[ARB]	# i: command string
#--
char	blank
int	ic, jc, kc, lc, it
pointer	sp, line

data	blank	/ ' ' /

begin
	call smark (sp)
	call salloc (line, maxch, TY_CHAR)

	jc = 0	# index of last char printed from command
	kc = 1	# one more than length of line
	lc = 0	# index of previous blank in line
	it = 1	# one more than number of times putfunc has been called

	for (ic = 1; command[ic] != EOS; ic = ic + 1) {
	    if (kc > maxch) {
		if (lc <= 1)
		    lc = kc

		Memc[line+lc-1] = EOS
		ic = jc + lc - 1
		jc = ic
		kc = 1
		lc = 0

		call putfunc (fd, it, Memc[line])
		it = it + 1

	    } else {
		if (command[ic] <= blank)
		    lc = kc

		Memc[line+kc-1] = command[ic]
		kc = kc + 1
	    }
	}

	if (kc > 1) {
	    Memc[line+kc-1] = EOS
	    call putfunc (fd, it, Memc[line])
	}

	call sfree (sp)
end

# PUTFILE -- Write the line to a text file

procedure putfile (fd, it, line)

int	fd		# i: file descriptor
int	it		# i: iteration of call (ignored)
char	line[ARB]	# i: line to add to output file
#--

begin
	call fprintf (fd, "%s\n")
	call pargstr (line)

end

# PUTHEAD -- Write the line to a table header

procedure puthead (tp, it, line)

pointer	tp		# i: table descriptor
int	it		# i: iteration of call
char	line[ARB]	# i: line to add to table header
#--
pointer	sp, name

begin
	call smark (sp)
	call salloc (name, SZ_FNAME, TY_CHAR)

	call makename ("expr", it, 2, Memc[name], SZ_FNAME)
	call tbhadt (tp, Memc[name], line)

	call sfree (sp)
end

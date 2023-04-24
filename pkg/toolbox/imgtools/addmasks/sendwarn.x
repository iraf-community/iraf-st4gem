# SEND_WARN -- Send a warning message about illegal mask value
#
# B.Simon	12-Apr-91	Original

procedure send_warn (index, iline, icol, value)

int	index		# i: Index mask name to print
int	iline		# i: Line number of bad pixel
int	icol		# i: Column number of mask pixel
int	value		# i: Value of mask pixel
pointer	file		# i: Input filename template
char	output[1]	# i: Output filename
#--
char	ofile[SZ_FNAME]
pointer	ifiles

int	junk
pointer sp, temp

string	badflag  "%s[%d,%d] has illegal value (%d). Treated as if 0.\n"

int	imtrgetim()

begin
	call smark (sp)
	call salloc (temp, SZ_FNAME, TY_CHAR)

	if (index == 0) {
	    call strcpy (ofile, Memc[temp], SZ_FNAME)
	} else {
	    junk = imtrgetim (ifiles, index, Memc[temp], SZ_FNAME)
	}

	call eprintf (badflag)
	call pargstr (Memc[temp])
	call pargi (icol)
	call pargi (iline)
	call pargi (value)

	call sfree (sp)
	return

# SET_WARN -- Set arguments of warning message

	entry set_warn (file, output)

	ifiles = file
	call strcpy (output, ofile, SZ_FNAME)

end
 

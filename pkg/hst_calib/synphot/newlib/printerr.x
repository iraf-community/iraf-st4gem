#* HISTORY *
#* B.Simon	28-Apr-94	original

# PRINTERR_INT -- Print an error message containing an integer value

procedure printerr_int (errmsg, ival)

char	errmsg[ARB]	# i: Error message
int	ival		# i: Integer value
#--
pointer	sp, line

begin
	call smark (sp)
	call salloc (line, SZ_LINE, TY_CHAR)

	call sprintf (Memc[line], SZ_LINE, "%s (%d)")
	call pargstr (errmsg)
	call pargi (ival)

	call error (1, Memc[line])
	call sfree (sp)

end

# PRINTERR_REAL -- Print an error message containing a real value

procedure printerr_real (errmsg, rval)

char	errmsg[ARB]	# i: Error message
real	rval		# i: Real value
#--
pointer	sp, line

begin
	call smark (sp)
	call salloc (line, SZ_LINE, TY_CHAR)

	call sprintf (Memc[line], SZ_LINE, "%s (%g)")
	call pargstr (errmsg)
	call pargr (rval)

	call error (1, Memc[line])
	call sfree (sp)

end

# PRINTERR_STR -- Print an error message containing a string

procedure printerr_str (errmsg, strval)

char	errmsg[ARB]	# i: Error message
char	strval[ARB]	# i: String text
#--
pointer	sp, line

begin
	call smark (sp)
	call salloc (line, SZ_LINE, TY_CHAR)

	call sprintf (Memc[line], SZ_LINE, "%s (%s)")
	call pargstr (errmsg)
	call pargstr (strval)

	call error (1, Memc[line])
	call sfree (sp)

end


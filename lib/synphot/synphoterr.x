include	"libsynphot.h"

# SYNPHOTERR -- Write a message and take an error exit from synphot

procedure synphoterr (errmsg, str)

char	errmsg[ARB]	# i: Error message
char	str[ARB]	# i: String to be appended to message
#--
errchk	synerrflag

begin
	call synerrflag (MISC_ERR, errmsg, str)
end

# SYNINTERR -- Write an error message with an integer value

procedure syninterr (errmsg, value)

char	errmsg[ARB]	# i: Error message
int	value		# i: integer value
#--
pointer	sp, str
errchk	synerrflag

begin
	call smark (sp)
	call salloc (str, SZ_FNAME, TY_CHAR)

	call sprintf (Memc[str], SZ_FNAME, "%d")
	call pargi (value)

	call synerrflag (MISC_ERR, errmsg, Memc[str])

	call sfree (sp)
end

# SYNREALERR -- Write an error message with an real value

procedure synrealerr (errmsg, value)

char	errmsg[ARB]	# i: Error message
real	value		# i: real value
#--
pointer	sp, str
errchk	synerrflag

begin
	call smark (sp)
	call salloc (str, SZ_FNAME, TY_CHAR)

	call sprintf (Memc[str], SZ_FNAME, "%f")
	call pargr (value)

	call synerrflag (MISC_ERR, errmsg, Memc[str])

	call sfree (sp)
end

# SYNERRFLAG -- A version that also allows you to set the error code

procedure synerrflag (code, errmsg, str)

int	code		# i: Error code 
char	errmsg[ARB]	# i: Error message
char	str[ARB]	# i: String to be appended to message
#--
char	line[SZ_LINE]

begin
	call strcpy (errmsg, line, SZ_LINE)
	call strcat (": ", line, SZ_LINE)
	call strcat (str, line, SZ_LINE)

	call error (code, line)
end


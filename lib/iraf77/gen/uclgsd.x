include	<iraf77.h>

# UCLGS? -- Get a scalar from the CL.
# If the parameter is list structured, the current list element is returned.

procedure uclgsd (f77par, buf, istat)

% 	character*(*) f77par
double	buf
int	istat

double	clgetd()
int	junk, clglpd()

char 	parnam[SZ_FNAME]
char 	parstr[SZ_FNAME]
char	partyp[SZ_FNAME]
int	stridx()
errchk	clgstr

begin
	istat = ER_OK
	# Convert character string to SPP string
	call f77upk (f77par, parnam, SZ_FNAME)
	# Get parameter's type from the cl (scalar, list, array)
	call strcpy (parnam, parstr, SZ_FNAME)
	call strcat (".p_xtype", parstr, SZ_FNAME)
	iferr (	call clgstr (parstr, partyp, SZ_FNAME) )  {
	      istat = ER_CLNOTFND
	      return
	}
	if (partyp[1] == '*') {
	   junk = clglpd (parnam, buf)	# list parameter
	   if (junk == EOF)
		istat = ER_CLEOF		# end of list
	   if (IS_INDEFD (buf))
		istat = ER_CLUNDF
	} else if (stridx ("a", partyp) > 0) {
	   istat = ER_CLBADTYP
	} else  {
	   buf = clgetd(parnam)		# simple scalar parameter
	   if (IS_INDEFD (buf))
		istat = ER_CLUNDF
	}
	return

end

include	<iraf77.h>

# UCLGST -- Get a character string from the CL.
# If the parameter is list structured, the current list element is returned.

procedure uclgst (f77par, f77buf, istat)

% 	character*(*) f77par
% 	character*(*) f77buf
int	istat

char	buf[SZ_LINE]
int	junk, clglstr()

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
	   junk = clglstr (parnam, buf, SZ_LINE)	# list parameter
	   if (junk == EOF)
		istat = ER_CLEOF		# end of list
	} else if (stridx ("a", partyp) > 0) {
	   istat = ER_CLBADTYP
	} else  {
	   call clgstr (parnam, buf, SZ_LINE)		# simple string
	}
	# Be sure the user gets a blank in case of an error
%       f77buf = ' '
	if (istat == ER_OK)
	   # convert string back to FORTRAN
	   call f77pak (buf, f77buf, SZ_LINE)
	return

end

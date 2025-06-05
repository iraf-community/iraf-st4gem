include	<iraf77.h>

# UCLPST -- Put a character string into the CL.

procedure uclpst (f77par, f77buf, istat)

% 	character*(*) f77par
% 	character*(*) f77buf
int	istat

char 	parnam[SZ_FNAME]
char 	bufnam[SZ_LINE]

begin
	istat = ER_OK
	# Convert parameter's name to SPP name
	call f77upk (f77par, parnam, SZ_FNAME)
	# Convert FORTRAN character string to SPP string
	call f77upk (f77buf, bufnam, SZ_LINE)
	# put the parameter
	iferr ( call clpstr (parnam, bufnam) )
	   istat = ER_CLPUT
	return

end

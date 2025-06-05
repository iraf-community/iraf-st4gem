include	<iraf77.h>

# UCLPS? -- Put a scalar into the CL.

procedure uclpsx (f77par, buf, istat)

% 	character*(*) f77par
complex	buf
int	istat

char 	parnam[SZ_FNAME]
# PIXEL   clput$t()

begin
	istat = ER_OK
	# Convert character string to SPP string
	call f77upk (f77par, parnam, SZ_FNAME)
	# put the parameter
	iferr ( call clputx (parnam, buf) )
	   istat = ER_CLPUT
	return

end

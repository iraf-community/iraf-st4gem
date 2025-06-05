include <iraf77.h>

# UCLRCU-- Return the next cursor value from a list structured cursor type
# parameter.

procedure uclrcu (f77par, wx, wy, key, f77str, istat)

% 	character*(*) f77par
% 	character*(*) f77str
real	wx, wy			# cursor coordinates
int	key			# keystroke value of cursor event
int	istat

int	cstat			# cursor return status
int	wcs			# wcs to which coordinates belong

pointer	sp, parnam, strval
int	txtlen			# size of input string

int	clgcur()

begin
	istat = ER_OK

	call smark (sp)

	# Allocate space for parameter name
	txtlen = len (f77par)
	call salloc (parnam, txtlen, TY_CHAR)

	# Convert parameter name to SPP char
	call f77upk (f77par, Memc[parnam], txtlen)

	# Allocate space for returned text
	txtlen = len (f77str)
	call salloc (strval, txtlen, TY_CHAR)

	# Initialize string value
	Memc[strval] = EOS

	# Read cursor
	cstat = clgcur (Memc[parnam], 
	    wx, wy, wcs, key, Memc[strval], txtlen)

	if (cstat == EOF)  {
	    # Reached the end of cursor list
	    istat = ER_CLEOF
	    return
	}

	if (Memc[strval] != EOS)
	    # Convert sting to FORTRAN CHARACTER
	    call f77pak (Memc[strval], f77str, SZ_LINE)
end

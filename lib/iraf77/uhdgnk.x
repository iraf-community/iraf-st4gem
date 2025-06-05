include <iraf77.h>

# UHDGNK -- Get the next field name matching the given template
# from an image header database. Sorting of the field list is
# optional. A prior call to UHDOKL is necessary to open the sorted
# list.

procedure uhdgnk (fn, f77out, istat)

pointer fn			# field name list descriptor
%	character*(*) f77out
int	istat			# error code indicator
int	nchar, imgnfn()
char	outpar[SZ_KEYWORD]

begin
	istat = ER_OK
	nchar = imgnfn (fn, outpar, SZ_KEYWORD)	
	if (nchar == EOF) 
	   istat = ER_EOF
	call f77pak (outpar, f77out, SZ_KEYWORD)
	return
end
	

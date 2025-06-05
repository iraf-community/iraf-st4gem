include <iraf77.h>

# UHDOKL -- open a sorted/unsorted field name list.

procedure uhdokl (im, f77tp, sortk, listp, istat)

pointer	im	# image descriptor
%	character*(*) f77tp
bool	sortk
pointer listp
int	istat, imofnl()
char	template[SZ_FNAME]

begin
	istat = ER_OK
	call f77upk (f77tp, template, SZ_FNAME)
	iferr ( listp = imofnl (im, template, sortk))
		istat = ER_BADINPAT
	return
end

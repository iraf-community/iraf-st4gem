include <iraf77.h>
# UIMRNM -- Rename an image.

procedure uimrnm (old, new, ier)

%	character*(*) old
%	character*(*) new
int	ier

pointer oldp, newp, sp

begin
	ier = ER_OK
	call smark (sp)
	call salloc (oldp, SZ_PATHNAME, TY_CHAR)
	call salloc (newp, SZ_PATHNAME, TY_CHAR)

	call f77upk (old, Memc[oldp], SZ_PATHNAME)
	call f77upk (new, Memc[newp], SZ_PATHNAME)
	iferr {
	   call imrename (Memc[oldp], Memc[newp])
	} then
	  ier = ER_IMRENAME

	call sfree (sp)
end

include <iraf77.h>
# UIMDEL -- Delete an image.

procedure uimdel (image, ier)

%	character*(*) image
int	ier

pointer sp, imp

begin
	ier = ER_OK
	call smark (sp)
	call salloc (imp, SZ_PATHNAME, TY_CHAR)

	call f77upk (image, Memc[imp], SZ_PATHNAME)
	iferr {
	   call imdelete (Memc[imp])
	} then 
	   ier = ER_IMDELETE

	call sfree (sp)
end

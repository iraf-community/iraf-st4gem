include <syserr.h>
include <imio.h>
include	<iraf77.h>

# UHDDSP -- Delete the named header keyword.

procedure uhddsp (im, keyw, ier)

pointer	im			# imfort image descriptor
%       character*(*) keyw
int	ier

pointer	sp, kp
int	errcode()

begin
	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	iferr (call imdelf (im, Memc[kp])) {
	    ier = errcode()
	    if (ier == SYS_IDBNODEL)
		ier = ER_HDNODEL
	    else
		ier = ER_HDDELNXKW
	} else {
	    ier = ER_OK
	    IM_UPDATE(im) = YES
	}

	call sfree (sp)
end

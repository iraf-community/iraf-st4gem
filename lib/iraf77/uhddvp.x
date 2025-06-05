include <syserr.h>
include <imio.h>
include	<iraf77.h>

# UHDDVP -- Delete the named array header keyword.

procedure uhddvp (im, keyw, felem, nelem, ier)

pointer	im			# imfort image descriptor
%       character*(*) keyw
int 	felem			# first element to delete
int 	nelem			# number of elements to delete
int	ier

pointer	sp, kp
int	errcode(), strlen(), itoc()
int	len, i, nch

begin
	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	len = strlen (Memc[kp])
	iferr {
	    for (i = felem; i < felem+nelem; i=i+1) {
	       nch = itoc (i, Memc[kp+len], SZ_KEYWORD)
	       call imdelf (im, Memc[kp])
	    }
	} then {
	    ier = errcode()
	    if (ier == SYS_IDBNODEL)
		ier = ER_HDNODEL
	    else
		ier = ER_OK	#No problem if key does not exist.
	} else {
	    ier = ER_OK
	    IM_UPDATE(im) = YES
	}

	call sfree (sp)
end

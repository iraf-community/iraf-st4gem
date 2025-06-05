include <syserr.h>
include	<iraf77.h>

# UHDGST-- Return the value of the named header keyword as a character
# string.

procedure uhdgst (im, keyw, sval, ier)

pointer	im			# imfort image descriptor
%       character*(*) keyw
%	character*(*) sval
int	ier
int     errcode()

pointer	sp, kp, vp

begin
	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)
	call salloc (vp, SZ_VALSTR, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call strlwr (Memc[kp])
	iferr (call imgstr (im, Memc[kp], Memc[vp], SZ_VALSTR)) {
	    ier = errcode()
	    if (ier == SYS_IDBKEYNF) 
	       ier = ER_HDRPARNF
	    else
	       ier = ER_HDRPARTY
	} else {
	    call f77pak (Memc[vp], sval, ARB)
	    ier = ER_OK
	}

	call sfree (sp)
end

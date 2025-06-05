include <syserr.h>
include <imio.h>
include	<iraf77.h>

# UHDPST -- Update the value of a keyword of type string.

procedure uhdpst (im, keyw, sval, ier)

pointer	im			# image descriptor
%       character*(*) keyw
%	character*(*) sval
int	ier

pointer	sp, kp, vp, cp
int	errcode(), idb_kwlookup()

begin
	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)
	call salloc (vp, SZ_VALSTR, TY_CHAR)
	call salloc (cp, SZ_VALSTR, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call f77upk (sval, Memc[vp], SZ_VALSTR)
	call strlwr (Memc[kp])
	if (idb_kwlookup (Memc[kp]) == 0)
	   call strupr (Memc[kp])
	iferr (call  impstr (im, Memc[kp], Memc[vp])) {
	    ier = errcode()
	    if (ier == SYS_IDBOVFL)
		ier = ER_HDRNOSP
	    else if (ier == SYS_IDBKEYNF)
	        ier = ER_HDRPARNF
	    else
	        ier = ER_HDRPARTY
	} else {
	    ier = ER_OK
	    IM_UPDATE(im) = YES
	}

	call sfree (sp)
end

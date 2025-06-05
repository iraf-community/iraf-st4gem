include <syserr.h>
include <imio.h>
include	<iraf77.h>

# UHDPS? -- Update a keyword value of the specified type.

procedure uhdpsi (im, keyw, ival, ier)

pointer	im			# image descriptor
%       character*(*) keyw
int	ival
int	ier

pointer	sp, kp, cp
int	errcode(), idb_kwlookup()

begin
	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)
	call salloc (cp, SZ_VALSTR, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call strlwr (Memc[kp])
	if (idb_kwlookup (Memc[kp]) == 0)
	   call strupr (Memc[kp])
	iferr (call imputi (im, Memc[kp], ival)) {
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

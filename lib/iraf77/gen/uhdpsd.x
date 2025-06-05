include <syserr.h>
include <imio.h>
include	<iraf77.h>

# UHDPS? -- Update a keyword value of the specified type.

procedure uhdpsd (im, keyw, dval, ier)

pointer	im			# image descriptor
%       character*(*) keyw
double	dval
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
	iferr (call imputd (im, Memc[kp], dval)) {
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

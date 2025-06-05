include <syserr.h>
include <imio.h>
include	<iraf77.h>

# UHDPV? -- Update a keyword value array of the specified type.

procedure uhdpvt (im, keyw, felem, nelem, tval, ier)

pointer	im			# image descriptor
%       character*(*) keyw
% 	character*(*) tval(1)
int	felem			# first vector element to update
int	nelem			# number of elements to update
int	ier

pointer	sp, kp, vp
int	errcode(), strlen(), itoc()
int	i, len, nch, idb_kwlookup()

begin
	if (felem <= 0) {
	   ier = ER_HDBADFELEM
	   return
	}
	if (nelem <= 0) {
	   ier = ER_HDBADNELEM
	   return
	}
	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)
	call salloc (vp, SZ_VALSTR, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call strlwr (Memc[kp])
	if (idb_kwlookup (Memc[kp]) == 0)
	   call strupr (Memc[kp])
        len = strlen (Memc[kp])
	iferr {
	    for (i = felem; i < felem+nelem; i=i+1) {
	       nch = itoc (i, Memc[kp+len], SZ_KEYWORD)
	       call f77upk (tval[i], Memc[vp], SZ_VALSTR)
	       call  impstr (im, Memc[kp], Memc[vp])	  
	    }
	} then {
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

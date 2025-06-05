include <syserr.h>
include <imio.h>
include	<iraf77.h>

# UHDPV? -- Update a keyword value array of the specified type.

procedure uhdpvl (im, keyw, felem, nelem, lval, ier)

pointer	im			# image descriptor
%       character*(*) keyw
long	lval[ARB]		# user array with new values
int	felem			# first vector element to update
int	nelem			# number of elements to update
int	ier

pointer	sp, kp
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

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call strlwr (Memc[kp])
	if (idb_kwlookup (Memc[kp]) == 0)
	   call strupr (Memc[kp])
	len = strlen (Memc[kp])
	iferr {
	    for (i = felem; i < felem+nelem; i=i+1) {
	       nch = itoc (i, Memc[kp+len], SZ_KEYWORD)
	       call imputl (im, Memc[kp], lval[i])
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

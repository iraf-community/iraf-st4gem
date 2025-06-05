include <syserr.h>
include	<iraf77.h>

# UHDGS? -- Return the value of the named header keyword 
#           in the specified type.

procedure uhdgss (im, keyw, sval, ier)

pointer	im			# imfort image descriptor
%       character*(*) keyw
short	sval
int	ier

pointer	sp, kp
short	imgets()
int	errcode()

begin
	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call strlwr (Memc[kp])
	iferr (sval = imgets (im, Memc[kp])) {
	    ier = errcode()
	    if (ier == SYS_IDBKEYNF)
	       ier = ER_HDRPARNF
	    else
	       ier = ER_HDRPARTY 
	} else
	    ier = ER_OK

	call sfree (sp)
end

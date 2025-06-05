include <syserr.h>
include	<iraf77.h>

# UHDGS? -- Return the value of the named header keyword 
#           in the specified type.

procedure uhdgsi (im, keyw, ival, ier)

pointer	im			# imfort image descriptor
%       character*(*) keyw
int	ival
int	ier

pointer	sp, kp
int	imgeti()
int	errcode()

begin
	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call strlwr (Memc[kp])
	iferr (ival = imgeti (im, Memc[kp])) {
	    ier = errcode()
	    if (ier == SYS_IDBKEYNF)
	       ier = ER_HDRPARNF
	    else
	       ier = ER_HDRPARTY 
	} else
	    ier = ER_OK

	call sfree (sp)
end

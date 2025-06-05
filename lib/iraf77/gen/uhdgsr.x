include <syserr.h>
include	<iraf77.h>

# UHDGS? -- Return the value of the named header keyword 
#           in the specified type.

procedure uhdgsr (im, keyw, rval, ier)

pointer	im			# imfort image descriptor
%       character*(*) keyw
real	rval
int	ier

pointer	sp, kp
real	imgetr()
int	errcode()

begin
	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call strlwr (Memc[kp])
	iferr (rval = imgetr (im, Memc[kp])) {
	    ier = errcode()
	    if (ier == SYS_IDBKEYNF)
	       ier = ER_HDRPARNF
	    else
	       ier = ER_HDRPARTY 
	} else
	    ier = ER_OK

	call sfree (sp)
end

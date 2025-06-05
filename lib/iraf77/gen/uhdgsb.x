include <syserr.h>
include	<iraf77.h>

# UHDGS? -- Return the value of the named header keyword 
#           in the specified type.

procedure uhdgsb (im, keyw, bval, ier)

pointer	im			# imfort image descriptor
%       character*(*) keyw
bool	bval
int	ier

pointer	sp, kp
bool	imgetb()
int	errcode()

begin
	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call strlwr (Memc[kp])
	iferr (bval = imgetb (im, Memc[kp])) {
	    ier = errcode()
	    if (ier == SYS_IDBKEYNF)
	       ier = ER_HDRPARNF
	    else
	       ier = ER_HDRPARTY 
	} else
	    ier = ER_OK

	call sfree (sp)
end

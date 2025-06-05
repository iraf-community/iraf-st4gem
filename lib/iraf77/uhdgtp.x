include	<iraf77.h>

# UHDGTP -- Get the datatype for a keyword.

procedure uhdgtp (im, keyw, dtype, ier)

pointer	im			# image descriptor
%	character*(*) keyw
int	dtype			# receives datatype code
int	ier

pointer	sp, kp, cp
int	imgftype()

begin
	ier = ER_OK

	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)
	call salloc (cp, SZ_VALSTR, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	iferr (dtype = imgftype (im, Memc[kp]))
	    ier = ER_HDRPARNF
	else
	    switch (dtype) {
	    # Map data type returned by imgftype into F77/VOS.  CB, 2/1/88
	    case TY_BOOL:
		dtype = TY_BOOL
	    case TY_CHAR:
		dtype = -1
	    case TY_INT:
		dtype = TY_LONG
	    case TY_LONG:
		dtype = TY_LONG		# TY_INT perhaps ?
	    case TY_REAL:
		dtype = TY_REAL
	    }

	call sfree (sp)
end

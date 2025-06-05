include <imio.h>
include	<iraf77.h>

# UHDAHS -- Add a history, comment, or blank (eight spaces) keyword string
# to a data file header.

procedure uhdahs (im, keyw, comm, ier)

pointer	im			# image descriptor
%       character*(*) keyw
%	character*(*) comm

int	ier

pointer	sp, kp, cp
errchk  imputh

begin
	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)
	call salloc (cp, SZ_CMDLINE-1, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call f77upk (comm, Memc[cp], SZ_CMDLINE-1)

	# note that the case of keyw does not matter; 
	# imputh will force upper case anyway

	iferr (call imputh (im, Memc[kp], Memc[cp])) 
	    ier = ER_HDRNOSP
	else {
 	    ier = ER_OK
	    IM_UPDATE(im) = YES
	}

	call sfree (sp)
end

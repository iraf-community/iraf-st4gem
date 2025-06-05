include <syserr.h>
include <imhdr.h>
include <imio.h>
include	<iraf77.h>

# UHDAS? -- Add a new keyword of the specified type. A comment string
# can be associated with the keyword. It is not an error if the parameter 
# already exists.

procedure uhdasd (im, keyw, dval, comm, htype, ier)

pointer	im			# imfort image descriptor
%       character*(*) keyw
double	dval
%       character*(*) comm
int	ier
int	htype

pointer	sp, kp, cp
int	errcode(), imaccf(), lenv, idb_kwlookup()
errchk  imadcf, imputd

begin
	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)
	call salloc (cp, SZ_VALSTR, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call f77upk (comm, Memc[cp], SZ_VALSTR)
	call strlwr (Memc[kp])
	if (idb_kwlookup (Memc[kp]) == 0)
	   call strupr (Memc[kp])	
	lenv = 0			# length for string type only
	if (imaccf (im, Memc[kp]) == NO)
	   iferr (call imadcf (im, Memc[kp], lenv, TY_DOUBLE, Memc[cp], htype)) {
	       ier = errcode()
               if (ier == SYS_IDBOVFL)
	         ier = ER_HDRNOSP
	       else 
	         ier = ER_HDRPARRDF
	       return
	   } 
	iferr (call imputd (im, Memc[kp], dval)) 
	              ier = ER_HDRPARTY
	else {
	   ier = ER_OK
	   IM_UPDATE(im) = YES
	}

	call sfree (sp)
end

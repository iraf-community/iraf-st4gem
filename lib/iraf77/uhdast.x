
include <syserr.h>
include <imio.h>
include	<iraf77.h>

# UHDAST -- Add a new keyword of type string to a data file header.

procedure uhdast (im, keyw, sval, comm, htype, ier)

pointer	im			# image descriptor
%       character*(*) keyw
%	character*(*) sval
%	character*(*) comm
int     htype			# imheader parameter flag for global
				# or group specific.
int	ier, strlen()
int 	lenv, ,idb_kwlookup()

pointer	sp, kp, vp, cp
int	errcode(), imaccf()
errchk  imadcf, impstr

begin
	call smark (sp)
	call salloc (kp, SZ_KEYWORD, TY_CHAR)
	call salloc (vp, SZ_VALSTR, TY_CHAR)
	call salloc (cp, SZ_VALSTR, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call f77upk (sval, Memc[vp], SZ_VALSTR)
	call f77upk (comm, Memc[cp], SZ_VALSTR)
	call strlwr (Memc[kp])
	if (idb_kwlookup (Memc[kp]) == 0)
	   call strupr(Memc[kp])
	lenv = strlen (Memc[vp])
	if (imaccf (im, Memc[kp]) == NO)
	    iferr (call imadcf (im, Memc[kp], lenv, TY_CHAR, Memc[cp], htype)) {
	       ier = errcode()
               if (ier == SYS_IDBOVFL)
	          ier = ER_HDRNOSP
	       else 
	          ier = ER_HDRPARRDF
	       return	 
            } 
	iferr (call impstr (im, Memc[kp], Memc[vp])) 
	           ier = ER_HDRPARTY
	else {
	   ier = ER_OK
	   IM_UPDATE(im) = YES
	}

	call sfree (sp)
end

include <syserr.h>
include <imhdr.h>
include <imio.h>
include	<iraf77.h>

# UHDAV? -- Add a new keyword array of the specified type. A comment string
# can be associated with the keyword. It is not an error if the parameter 
# already exists.

procedure uhdavi (im, keyw, felem, nelem, ival, comm, htype, ier)

pointer	im			# image descriptor
%       character*(*) keyw
int	ival[ARB]
%       character*(*) comm
int     felem			# first element to add
int	nelem			# number of elements to add
int	htype			# add to header file (0) or to gpb (1)
int	ier

pointer	sp, kp, cp
int	errcode(), imaccf(), itoc(), strlen()
int	len, i, nch, lenv, idb_kwlookup()

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
	call salloc (cp, SZ_VALSTR, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call f77upk (comm, Memc[cp], SZ_VALSTR)
	call strlwr (Memc[kp])
	if (idb_kwlookup (Memc[kp]) == 0)
	   call strupr (Memc[kp])
	len = strlen (Memc[kp])
	lenv = 0                 # length for string type only
	iferr {
	   for (i = felem; i < felem+nelem; i=i+1) {
	      nch = itoc (i, Memc[kp+len], SZ_KEYWORD)
	      if (imaccf (im, Memc[kp]) == NO) {
	          call imadcf (im, Memc[kp], lenv, TY_INT, Memc[cp], htype)
	      } else {
	          ier = ER_HDRPARRDF
	          goto 99
	      }
	      call imputi (im, Memc[kp], ival[i])
	    }
	} then {
	    ier = errcode()
	    if (ier == SYS_IDBOVFL)
	       ier = ER_HDRNOSP
	    else if (ier == SYS_IDBREDEF)
	       ier = ER_HDRPARRDF
	    else
	       ier = ER_HDRPARTY
	} else {
	    ier = ER_OK
	    IM_UPDATE(im) = YES	
	}
99
	call sfree (sp)
end

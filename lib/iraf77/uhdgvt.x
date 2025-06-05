include <syserr.h>
include	<iraf77.h>

# UHDGVT -- Return the values of the named header keyword 
#           in the string type vector.

procedure uhdgvt (im, keyw, felem, nelem, nvals, tval, ier)

pointer	im			# imfort image descriptor
%       character*(*) keyw
%	character*(*) tval(1)
int	felem			# index of first vector to get
int	nelem			# number of elements to get
int	nvals			# number of elements retrieved
int	ier

pointer	sp, kp, vp
int	errcode(), strlen(), itoc()
int     i, len, nch
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
	call salloc (vp, SZ_VALSTR, TY_CHAR)

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call strlwr (Memc[kp])

	# get first element 
	nvals = 0
	for (i = felem; i < felem+nelem; i=i+1) {
	    if ( i == felem )
	       len = strlen (Memc[kp])
	    nch = itoc (i, Memc[kp+len], SZ_KEYWORD)
	    iferr (call imgstr (im, Memc[kp], Memc[vp], SZ_VALSTR)) {
	       ier = errcode()
	       if (ier == SYS_IDBKEYNF)
		  if (nvals == 0) 
		     ier = ER_HDRPARNF
		  else
		     ier = ER_OK
	       else
	          ier = ER_HDRPARTY 
	       goto 99
	    } else 
	        call f77pak (Memc[vp], tval[i], SZ_VALSTR)
	        nvals = nvals + 1
	}
        ier = ER_OK
99
	call sfree (sp)
end

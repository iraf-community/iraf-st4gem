include <syserr.h>
include <iraf77.h>

# UHDGV? -- Return the values of the named header keyword 
#           in the specific data type vector.

procedure uhdgvd (im, keyw, felem, nelem, nvals, dval, ier)

pointer	im			# imfort image descriptor
%       character*(*) keyw
double	dval[ARB]
int	felem			# index of first vector to get
int	nelem			# number of elements to get
int	nvals			# number of elements retrieved
int	ier

pointer	sp, kp
double	imgetd()
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

	call f77upk (keyw, Memc[kp], SZ_KEYWORD)
	call strlwr (Memc[kp])

	# get first element 
	nvals = 0
	for (i = felem; i < felem+nelem; i=i+1) {
	    if ( i == felem )
	       len = strlen (Memc[kp])
	    nch = itoc (i, Memc[kp+len], SZ_KEYWORD)
	    iferr (dval[i] = imgetd (im, Memc[kp])) {
	       ier = errcode()
	       if (ier == SYS_IDBKEYNF)
		  if (nvals == 0) 
		     ier = ER_HDRPARNF
		  else
		     ier = ER_OK
	       else
	          ier = ER_HDRPARTY 
	       goto 99
            }
	    nvals = nvals + 1
	}
        ier = ER_OK
99
	call sfree (sp)
end

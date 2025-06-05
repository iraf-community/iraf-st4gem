include	<iraf77.h>
include	<syserr.h>

# UIMOPN -- Opens an existing image (most probably from an SDAS group file)
# f77nam is the name returned from a call to the template expander and should
# have the form: filename[n][section] where n is the group to access.

procedure uimopn (f77nam, acc_mod, im_id, istat)

% 	character*(*) f77nam
int	acc_mod			# access mode (RO, RW)
pointer im_id			# image descriptor
int	istat

char 	imname[SZ_PATHNAME]
pointer immap()
int	ierr
int	errcode()

begin

	if ( acc_mod != READ_ONLY  && acc_mod != READ_WRITE)  {
	   istat = ER_IMBADACCMOD
	   return
	}
	istat = ER_OK
	# Convert character string to SPP string
	call f77upk (f77nam, imname, SZ_PATHNAME)
# The next statement is commented to allow for image names with 
# upper case. Nelson Zarate March 24 1993
#	call strlwr (imname)

	iferr ( im_id = immap (imname, acc_mod, 0) )  {
	   ierr = errcode ()
	   if (ierr == SYS_IMMAGNCPY)
		istat = ER_IMNOTIMAG
	   if (ierr == SYS_IMSYNSEC || ierr == SYS_IMDIMSEC ||
	       ierr == SYS_IMSTEPSEC )
		istat = ER_IMBADSEC
	   if (ierr == SYS_IKIEXTN)
	        istat = ER_IMBADEXTN
	   if (ierr == SYS_IKIOPEN)
  	       istat = ER_IMOPOLD
	   # reset pointer to zero to avoid passing back an invalid pointer
	   im_id = 0
	}
	return

end

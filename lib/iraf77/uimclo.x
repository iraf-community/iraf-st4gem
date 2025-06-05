include	<iraf77.h>

# UIMCLO -- Close all files related to a given image.

procedure uimclo (im_id, istat)

pointer im_id			# image descriptor
int	istat			# iraf77 error code
errchk	imunmap

begin

	istat = ER_OK

	iferr ( call imunmap (im_id) )
	     istat = ER_IMCLOS
	return

end

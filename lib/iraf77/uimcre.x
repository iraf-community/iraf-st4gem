define  DEF_USERAREA    10000	# Size of user area (in Header)
include	<imhdr.h>
include	<iraf77.h>
include <syserr.h>

# UIMCRE -- Create an image from scratch.
# f77nam is the name of the new image and is of the form:
# filename[p/n][section] where n is the number of images in the group-file.

procedure uimcre (f77nam, dtype, naxis, dimen, im_id, istat)

%	character*(*) f77nam
int	dtype			# image data_type
int	naxis			# image dimensionnality
long	dimen[IM_MAXDIM]	# size in each dimension
pointer im_id
int	istat

char 	imname[SZ_PATHNAME]
pointer immap()
int	i, ierr
int	errcode()

begin
	istat = ER_OK
	# Verify image size and datatype operands.
	if (naxis < 1 || naxis > IM_MAXDIM)
	    istat = ER_IMBADNAXIS
	if (istat == ER_OK)
	    do i = 1, naxis
		if (dimen[i] < 1)
		    istat = ER_IMBADDIMEN
	if (istat == ER_OK)
	    if (dtype != TY_USHORT && dtype != TY_SHORT && dtype != TY_REAL
		&& dtype != TY_INT && dtype != TY_LONG  && dtype != TY_DOUBLE
		&& dtype != TY_COMPLEX)
		   istat = ER_IMBADTYPE
	if (istat != ER_OK) 
	    return

	# Convert character string to SPP string
	call f77upk (f77nam, imname, SZ_PATHNAME)

	iferr ( im_id = immap (imname, NEW_IMAGE, DEF_USERAREA)) {
	   ierr = errcode ()
	   if (ierr == SYS_IMSECTNEWIM)
		istat = ER_IMILLSEC
	   if (ierr == SYS_IKIEXTN)
	        istat = ER_IMBADEXTN
	   if (ierr == SYS_IKIOPEN)
  	       istat = ER_IMOPNEW
	   # reset pointer to zero to avoid passing back an invalid pointer
	   im_id = 0
	   return
	}
	# Set pixel file characteristics
	call amovl (dimen, IM_LEN(im_id,1), naxis)
	IM_NDIM(im_id) = naxis
	IM_PIXTYPE(im_id) = dtype
	return

end

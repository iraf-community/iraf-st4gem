include	<imhdr.h>
include	<iraf77.h>

# UIMGMP -- Gets the data portion of the selected image into virtual memory
# and return an offset into a buffer for the selected data type.

procedure uimgmp (im_id, offset, istat)

pointer im_id			# image descriptor
pointer offset			# pointer into virtual memory
int	istat			# return status

int	naxis			# image dimensionnality
long	vs[IM_MAXDIM]		# section starting vector
long	ve[IM_MAXDIM]		# section ending vector
int	i
pointer sp, msg


pointer imggss()

pointer imggsi()

pointer imggsl()

pointer imggsr()

pointer imggsd()

pointer imggsx()


begin

	istat = ER_OK
	# Set up section to get the whole image
	naxis = IM_NDIM(im_id)
	do i = 1, naxis  {
	   vs[i] = 1
	   ve[i] = IM_LEN(im_id,i)
	}

	switch (IM_PIXTYPE(im_id)) {

	case TY_SHORT:
	# get the image into virtual memory.
	iferr ( offset = imggss (im_id, vs, ve, naxis) )
	   istat = ER_IMOFFSET

	case TY_INT:
	# get the image into virtual memory.
	iferr ( offset = imggsi (im_id, vs, ve, naxis) )
	   istat = ER_IMOFFSET

	case TY_LONG:
	# get the image into virtual memory.
	iferr ( offset = imggsl (im_id, vs, ve, naxis) )
	   istat = ER_IMOFFSET

	case TY_REAL:
	# get the image into virtual memory.
	iferr ( offset = imggsr (im_id, vs, ve, naxis) )
	   istat = ER_IMOFFSET

	case TY_DOUBLE:
	# get the image into virtual memory.
	iferr ( offset = imggsd (im_id, vs, ve, naxis) )
	   istat = ER_IMOFFSET

	case TY_COMPLEX:
	# get the image into virtual memory.
	iferr ( offset = imggsx (im_id, vs, ve, naxis) )
	   istat = ER_IMOFFSET

	default:
	   # error handling for unsupported datatype
	   call smark (sp)
	   call salloc (msg, SZ_LINE, TY_CHAR)
	   call sprintf (Memc[msg], SZ_LINE,
			    "uimgmp - unsupported data type %d")
	   call pargi (IM_PIXTYPE(im_id))
	   call error (0, Memc[msg])
	   call sfree (sp)
	}
end

include	<imhdr.h>
include	<iraf77.h>

# UIMGID -- Returns basic file description information for an open image.

procedure uimgid (im_id, dtype, naxis, dimen, istat)

pointer im_id			# image descriptor
int	dtype			# image data_type
int	naxis			# image dimensionnality
long	dimen[IM_MAXDIM]	# size in each dimension
int	istat

begin

	# Return image characteristics
	call amovl (IM_LEN(im_id,1), dimen, IM_MAXDIM)
	naxis = IM_NDIM(im_id)
	dtype = IM_PIXTYPE(im_id)
	istat = ER_OK
	return

end

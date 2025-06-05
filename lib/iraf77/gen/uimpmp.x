include	<imhdr.h>
include	<iraf77.h>

# UIMPMP -- Puts the data portion of the selected (output) image into an IMIO
# buffer.   This buffer will eventually be flushed.
# Simulates the VAX %VAL construct.

procedure uimpmp (im_id, in_offset, istat)

pointer im_id			# image descriptor (of output image)
pointer in_offset		# pointer into virtual memory (input image)
int	istat			# return status

int	naxis			# image dimensionnality
int 	npix			# image size in pixels
long	vs[IM_MAXDIM]		# section starting vector
long	ve[IM_MAXDIM]		# section ending vector
int	i
pointer out_offset, sp, msg


pointer impgss()

pointer impgsi()

pointer impgsl()

pointer impgsr()

pointer impgsd()

pointer impgsx()


begin

	istat = ER_OK
	# Set up section to put the whole image
	# also compute number of pixels in image
	naxis = IM_NDIM(im_id)
	npix = 1
	do i = 1, naxis  {
	   vs[i] = 1
	   ve[i] = IM_LEN(im_id,i)
	   npix = npix * ve[i]
	}

	switch (IM_PIXTYPE(im_id)) {

	case TY_SHORT:
	# put the image (in fact get an IMIO buffer)
	iferr ( out_offset = impgss (im_id, vs, ve, naxis) )
	   istat = ER_IMOFFSET
	# transfer from input buffer pointed to by in_offset to output buffer
	# pointed to by out_offset
	call amovs (Mems[in_offset], Mems[out_offset], npix)

	case TY_INT:
	# put the image (in fact get an IMIO buffer)
	iferr ( out_offset = impgsi (im_id, vs, ve, naxis) )
	   istat = ER_IMOFFSET
	# transfer from input buffer pointed to by in_offset to output buffer
	# pointed to by out_offset
	call amovi (Memi[in_offset], Memi[out_offset], npix)

	case TY_LONG:
	# put the image (in fact get an IMIO buffer)
	iferr ( out_offset = impgsl (im_id, vs, ve, naxis) )
	   istat = ER_IMOFFSET
	# transfer from input buffer pointed to by in_offset to output buffer
	# pointed to by out_offset
	call amovl (Meml[in_offset], Meml[out_offset], npix)

	case TY_REAL:
	# put the image (in fact get an IMIO buffer)
	iferr ( out_offset = impgsr (im_id, vs, ve, naxis) )
	   istat = ER_IMOFFSET
	# transfer from input buffer pointed to by in_offset to output buffer
	# pointed to by out_offset
	call amovr (Memr[in_offset], Memr[out_offset], npix)

	case TY_DOUBLE:
	# put the image (in fact get an IMIO buffer)
	iferr ( out_offset = impgsd (im_id, vs, ve, naxis) )
	   istat = ER_IMOFFSET
	# transfer from input buffer pointed to by in_offset to output buffer
	# pointed to by out_offset
	call amovd (Memd[in_offset], Memd[out_offset], npix)

	case TY_COMPLEX:
	# put the image (in fact get an IMIO buffer)
	iferr ( out_offset = impgsx (im_id, vs, ve, naxis) )
	   istat = ER_IMOFFSET
	# transfer from input buffer pointed to by in_offset to output buffer
	# pointed to by out_offset
	call amovx (Memx[in_offset], Memx[out_offset], npix)

	default:
	   # error handling for unsupported datatype
	   call smark (sp)
	   call salloc (msg, SZ_LINE, TY_CHAR)
	   call sprintf (Memc[msg], SZ_LINE,
			    "uimpmp - unsupported data type %d")
	   call pargi (IM_PIXTYPE(im_id))
	   call error (0, Memc[msg])
	   call sfree (sp)
	}
	return

end

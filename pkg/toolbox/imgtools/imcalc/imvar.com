# IMVAR.COM -- Common block used to hold image variables

pointer	imptr		# Image descriptors
pointer	dataptr		# Pointer to pixel array
pointer	typeptr		# Pointer to type array
pointer	gindex		# Pointer to group index array
pointer	gcount		# Pointer to group count array
pointer	oldline		# Previous line in image
pointer	newline		# Current line in image
pointer	outroot		# Output image rootname
pointer	outsect		# Output image section

int	nimage		# Number of images
int	npix		# Number of pixels in a line
int	nline		# Number of lines read from image

common	/imvar/	imptr, dataptr, typeptr, gindex, gcount,
		oldline, newline, outroot, outsect,
		nimage, npix, nline

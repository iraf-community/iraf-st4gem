#* HISTORY *
#* B.Simon	05-Mar-93	original
#* B.Simon	25-Jul-97	added line to set pixtype

include	<imhdr.h>

# WRTONEDIM -- Wrte a one dimensional image

procedure wrtonedim (image, crval, cd, ctype, spec, len)

char	image[ARB]	# i: spectrum name
real	crval		# i: world coordinate of first pixel
real	cd		# i: first derivative of world coordinates
char	ctype[ARB]	# i: world coordinate units
real	spec[ARB]	# i: spectrum values
int	len		# i: spectrum length
#--
long	aeon
pointer	im, wbuf

data	aeon  / 0 /

long	clktime()
pointer	immap(), impl1r()

begin
	# Open image and set dimensions

	im = immap (image, NEW_FILE, 0)
	IM_NDIM(im) = 1
	IM_LEN(im,1) = len
	IM_PIXTYPE(im) = TY_REAL

	# Write spectrum to image

	wbuf = impl1r (im)
	call amovr (spec, Memr[wbuf], len)

	# Set limits and write world coordinates

	call alimr (Memr[wbuf], len, IM_MIN(im), IM_MAX(im))
	IM_LIMTIME(im) = clktime (aeon)

	call imaddr (im, "CRPIX1", 1.0)
	call imaddr (im, "CRVAL1", crval)
	call imastr (im, "CTYPE1", ctype)
	call imaddr (im, "CD1_1", cd)

	# Close image

	call imunmap (im)
end

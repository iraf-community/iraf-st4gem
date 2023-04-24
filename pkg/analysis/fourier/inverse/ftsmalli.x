include <imhdr.h>
include "../fourier.h"

# ft_small_i -- inverse transform for small images
# This routine gets pointers to input & output image data and calls
# a routine to do the inverse Fourier transform in-memory.

procedure ft_small_i (fti, fto, decenter)

pointer fti		# i: pointer to FT structure for input image
pointer fto		# i: pointer to FT structure for output image
bool	decenter	# i: decenter output?
#--
pointer ximg		# scratch for complex copy of image
pointer iimRe, iimIm	# real & imaginary input image header pointers
pointer oimRe, oimIm	# real & imaginary output image header pointers
complex xnpts		# number of pixels in array (for normalization)
real	rnpts		# number of pixels in array
int	npts		# number of pixels in array
int	npix1, npix2	# lengths of first & second axes
bool	fwd		# forward transform?  (no)
bool	center		# center the output?  (no)

begin
	iimRe = FT_REPT(fti)
	iimIm = FT_IMPT(fti)
	oimRe = FT_REPT(fto)
	oimIm = FT_IMPT(fto)

	npix1 = IM_LEN(FT_IMAGE(fti),1)
	if (IM_NDIM(FT_IMAGE(fti)) > 1)
	    npix2 = IM_LEN(FT_IMAGE(fti),2)
	else
	    npix2 = 1

	npts = npix1 * npix2
	rnpts = real (npts)
	xnpts = complex (rnpts, 0.0)		# for normalization

	# Allocate space for a complex array to contain the entire image.
	call malloc (ximg, npts, TY_COMPLEX)

	# Get the input data and convert to complex.  Decenter if appropriate.
	call ft_ix_copy (fti, Memx[ximg], npix1, npix2, decenter)

	# Do the inverse Fourier transform in-place in the complex array.
	fwd = false
	call ft_cmplx (Memx[ximg], npix1, npix2, fwd)

	# Normalize by dividing by the total number of pixels.
	call adivkx (Memx[ximg], xnpts, Memx[ximg], npts)

	# Convert from complex, and write the data to output.
	center = false
	call ft_xo_copy (fto, Memx[ximg], npix1, npix2, center)

	call mfree (ximg, TY_COMPLEX)
end

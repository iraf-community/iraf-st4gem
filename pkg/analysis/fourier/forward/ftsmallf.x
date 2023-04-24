include <imhdr.h>
include "../fourier.h"

# ft_small_f -- forward transform for small images
# This routine gets pointers to input & output image data and calls
# a routine to do the forward Fourier transform in-memory.
#
# Phil Hodge, 19-Apr-1992  Subroutine created.
# Phil Hodge, 23-Jul-1993  Remove section that modified FT_CRPIX.

procedure ft_small_f (fti, fto, center)

pointer fti		# i: pointer to FT structure for input image
pointer fto		# i: pointer to FT structure for output image
bool	center		# i: center output?
#--
pointer ximg		# scratch for complex copy of image
pointer iimRe, iimIm	# real & imaginary input image header pointers
pointer oimRe, oimIm	# real & imaginary output image header pointers
int	npix1, npix2	# lengths of first & second axes
bool	fwd		# forward transform?  (yes)
bool	decenter	# decenter the input?  (no)

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

	# Allocate space for a complex array to contain the entire image.
	call malloc (ximg, npix1 * npix2, TY_COMPLEX)

	# Get the input data and convert to complex.
	decenter = false
	call ft_ix_copy (fti, Memx[ximg], npix1, npix2, decenter)

	# Do the forward Fourier transform in-place in the complex array.
	fwd = true
	call ft_cmplx (Memx[ximg], npix1, npix2, fwd)

	# Convert from complex, and write the data to output.
	call ft_xo_copy (fto, Memx[ximg], npix1, npix2, center)

	call mfree (ximg, TY_COMPLEX)
end

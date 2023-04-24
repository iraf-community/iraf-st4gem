include <imhdr.h>
include "../fourier.h"

# ft_small_x -- cross correlation for small images
# This routine reads two images into memory, takes the Fourier transform
# of each, multiplies the complex conjugate of the first by the second,
# takes the inverse Fourier transform, and writes the inverse transform
# to an output image.
#
# Phil Hodge,  2-Jul-1993  Subroutine created.
# Phil Hodge, 20-Jul-1993  Take complex conjugate of the first array instead
#			of the second array.

procedure ft_small_x (fti1, fti2, fto, naxis1, naxis2, center)

pointer fti1		# i: pointer to FT structure for input image
pointer fti2		# i: pointer to FT structure for PSF image
pointer fto		# i: pointer to FT structure for output image
int	naxis1, naxis2	# i: size of array used for Fourier transform
bool	center		# i: center output?
#--
pointer sp
pointer ximg1		# scratch for complex copy of image
pointer ximg2		# scratch for complex copy of PSF image
complex	xnpix		# naxis1 * naxis2, converted to complex
bool	fwd		# forward transform?
bool	decenter	# decenter the input?  (no)
errchk	salloc, ft_ix_copy

begin
	call smark (sp)

	# Allocate space for complex arrays to contain the images.
	call salloc (ximg1, naxis1 * naxis2, TY_COMPLEX)
	call salloc (ximg2, naxis1 * naxis2, TY_COMPLEX)

	# Get the input data and convert to complex.
	decenter = false
	call ft_ix_copy (fti1, Memx[ximg1], naxis1, naxis2, decenter)
	call ft_ix_copy (fti2, Memx[ximg2], naxis1, naxis2, decenter)

	# Do the forward Fourier transform in-place in the complex arrays.
	fwd = true
	call ft_cmplx (Memx[ximg1], naxis1, naxis2, fwd)
	call ft_cmplx (Memx[ximg2], naxis1, naxis2, fwd)

	# Take the complex conjugate of the first array.
	call acjgx (Memx[ximg1], Memx[ximg1], naxis1*naxis2)

	# Multiply the complex arrays, writing back to ximg1.
	call amulx (Memx[ximg1], Memx[ximg2], Memx[ximg1], naxis1*naxis2)

	# Normalize, working in-place in ximg1.
	xnpix = complex (naxis1*naxis2, 0)	# both arguments are integer
	call adivkx (Memx[ximg1], xnpix, Memx[ximg1], naxis1*naxis2)

	# Take the inverse Fourier transform in-place in ximg1.
	fwd = false
	call ft_cmplx (Memx[ximg1], naxis1, naxis2, fwd)

	# Convert from complex, and write the data to output.
	call ft_xo_copy (fto, Memx[ximg1], naxis1, naxis2, center)

	call sfree (sp)
end

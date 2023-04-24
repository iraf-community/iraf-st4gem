# This file contains ft_ix_copy, ft_xo_copy and ft_shift_x.
# The first two copy from input images to a complex array and from
# a complex array to output images respectively.  The rows and columns
# may optionally be shifted during the copy process.  The ft_shift_x
# routine shifts a 1-D complex array.
#
#  *****
# Note that the decenter option in ft_ix_copy requires the image and
# array to be the same size.  If they are not the same size, then
# decenter must be false.
#  *****
#
# Phil Hodge, 19-Jun-1992  Subroutines created.
# Phil Hodge,  2-Jul-1993  Array size nx by ny can be larger than image.

include <imhdr.h>
include "../fourier.h"

# ft_ix_copy -- copy input to a complex array
# This routine reads data from the input images, the real and/or imaginary
# parts, and copies to a complex array.
#
# If decenter is true, the image will be shifted so that the reference
# pixel is moved to [1,1].  In this case, the output array ximg must be
# the SAME SIZE as the image.

procedure ft_ix_copy (fti, ximg, nx, ny, decenter)

pointer fti			# i: pointer to FT structure for input image
complex ximg[nx,ny]		# o: array to contain image data
int	nx, ny			# i: array size (may be larger than image)
bool	decenter		# i: decenter input?
#--
pointer sp
pointer xline			# scratch for one line of image(s)
pointer iimRe, iimIm		# real & imaginary input image header pointers
pointer iRe, iIm		# pointers to real & imaginary input data
complex zero			# zero
int	ndim			# dimension of image
int	npix1, npix2		# size of image
int	shift[2]		# amount to shift for decentering
int	i, j
int	out			# if decenter=true, move line j to line out
pointer imgl2r()

begin
	zero = (0.0, 0.0)

	iimRe = FT_REPT(fti)
	iimIm = FT_IMPT(fti)

	ndim = IM_NDIM(FT_IMAGE(fti))
	npix1 = IM_LEN(FT_IMAGE(fti),1)
	if (ndim > 1)
	    npix2 = IM_LEN(FT_IMAGE(fti),2)
	else
	    npix2 = 1

	shift[1] = 0				# initial values
	shift[2] = 0

	if (decenter) {
	    shift[1] = 1 - nint (FT_CRPIX(fti,1))
	    if (npix2 > 1)
		shift[2] = 1 - nint (FT_CRPIX(fti,2))
	}

	# We only need to use this section if the shift is non-zero.
	if (shift[1] != 0 || shift[2] != 0) {

	    if (nx != npix1 || ny != npix2)
		call error (1, "ft_ix_copy:  image and array differ in size")

	    call smark (sp)
	    # We will copy an image line into this array, converting to
	    # complex, and then decenter this array by shift[1] and put
	    # the output into the appropriate line of ximg.
	    call salloc (xline, npix1, TY_COMPLEX)

	    # We will move line 1 of the image to ximg[1,out].
	    if (ndim == 1) {
		out = 1
	    } else {
		out = 1 + shift[2]
		out = mod (out, npix2)
		if (out < 0)
		    out = out + npix2
	    }

	    # Get the input data and convert to complex.
	    if (FT_REAL(fti) == YES && FT_IMAG(fti) == YES) {
		do j = 1, npix2 {
		    iRe = imgl2r (iimRe, j)
		    iIm = imgl2r (iimIm, j)
		    call apkxr (Memr[iRe], Memr[iIm], Memx[xline], npix1)
		    call ft_shift_x (Memx[xline], shift[1], ximg[1,out], npix1)
		    out = out + 1
		    if (out > npix2)
			out = out - npix2
		}
	    } else if (FT_REAL(fti) == YES) {
		do j = 1, npix2 {
		    iRe = imgl2r (iimRe, j)
		    call achtrx (Memr[iRe], Memx[xline], npix1)
		    call ft_shift_x (Memx[xline], shift[1], ximg[1,out], npix1)
		    out = out + 1
		    if (out > npix2)
			out = out - npix2
		}
	    } else if (FT_IMAG(fti) == YES) {
		do j = 1, npix2 {
		    iIm = imgl2r (iimIm, j)
		    do i = 0, npix1-1			# zero indexed
			Memx[xline+i] = complex (0.0, Memr[iIm+i])
		    call ft_shift_x (Memx[xline], shift[1], ximg[1,out], npix1)
		    out = out + 1
		    if (out > npix2)
			out = out - npix2
		}
	    }
	    call sfree (sp)

	} else {				# just copy without decentering

	    # Get the input data and convert to complex.
	    if (FT_REAL(fti) == YES && FT_IMAG(fti) == YES) {
		do j = 1, npix2 {
		    iRe = imgl2r (iimRe, j)
		    iIm = imgl2r (iimIm, j)
		    call apkxr (Memr[iRe], Memr[iIm], ximg[1,j], npix1)
		    do i = npix1+1, nx
			ximg[i,j] = zero
		}
	    } else if (FT_REAL(fti) == YES) {
		do j = 1, npix2 {
		    iRe = imgl2r (iimRe, j)
		    call achtrx (Memr[iRe], ximg[1,j], npix1)
		    do i = npix1+1, nx
			ximg[i,j] = zero
		}
	    } else if (FT_IMAG(fti) == YES) {
		do j = 1, npix2 {
		    iIm = imgl2r (iimIm, j)
		    do i = 1, npix1
			ximg[i,j] = complex (0.0, Memr[iIm+i-1])
		    do i = npix1+1, nx
			ximg[i,j] = zero
		}
	    }
	    do j = npix2+1, ny
		call amovkx (zero, ximg[1,j], nx)
	}
end

# ft_xo_copy -- copy complex array to output
# This routine copies from a complex array to output images, the real and/or
# imaginary parts.
#
# If center is true, the array will be shifted so pixel [1,1] is moved to
# the center of the image.
#
# The image may be smaller than the array, even if center is true.

procedure ft_xo_copy (fto, ximg, nx, ny, center)

pointer fto			# i: pointer to FT structure for output image
complex ximg[nx,ny]		# i: array containing data
int	nx, ny			# i: array size (may be larger than image)
bool	center			# i: center the output?
#--
pointer sp
pointer xline			# scratch for one line of image(s)
pointer oimRe, oimIm		# real & imaginary output image header pointers
pointer oRe, oIm		# pointers to real & imaginary output data
int	ndim			# dimension of image
int	npix1, npix2		# size of image
int	shift[2]		# amount to shift for centering
int	i, j
int	in			# if center=true, move line in to line j
pointer impl2r()

begin
	oimRe = FT_REPT(fto)
	oimIm = FT_IMPT(fto)

	ndim = IM_NDIM(FT_IMAGE(fto))
	npix1 = IM_LEN(FT_IMAGE(fto),1)
	if (ndim > 1)
	    npix2 = IM_LEN(FT_IMAGE(fto),2)
	else
	    npix2 = 1

	# We use this section if output is to be centered.
	if (center) {

	    shift[1] = npix1 / 2
	    shift[2] = npix2 / 2

	    call smark (sp)
	    # We will copy a line into scratch, shifting by shift[1],
	    # and unpack the scratch array into the output image(s).
	    call salloc (xline, nx, TY_COMPLEX)

	    # We will move line ximg[1,in] to line j of the image.
	    # This initial value of 'in' corresponds to j = 1.
	    if (ndim == 1) {
		in = 1
	    } else {
		in = 1 - shift[2] + ny
		in = mod (in, ny)
	    }

	    # Convert from complex, and write the data to output.
	    if (FT_REAL(fto) == YES && FT_IMAG(fto) == YES) {
		# top and bottom parts of ximg
		do j = 1, npix2 {
		    # Note that we use nx for ft_shift_x but npix1 for aupxr.
		    call ft_shift_x (ximg[1,in], shift[1], Memx[xline], nx)
		    oRe = impl2r (oimRe, j)
		    oIm = impl2r (oimIm, j)
		    call aupxr (Memx[xline], Memr[oRe], Memr[oIm], npix1)
		    in = in + 1
		    if (in > ny)
			in = in - ny
		}
	    } else if (FT_REAL(fto) == YES) {
		do j = 1, npix2 {
		    call ft_shift_x (ximg[1,in], shift[1], Memx[xline], nx)
		    oRe = impl2r (oimRe, j)
		    call achtxr (Memx[xline], Memr[oRe], npix1)
		    in = in + 1
		    if (in > ny)
			in = in - ny
		}
	    } else if (FT_IMAG(fto) == YES) {
		do j = 1, npix2 {
		    call ft_shift_x (ximg[1,in], shift[1], Memx[xline], nx)
		    oIm = impl2r (oimIm, j)
		    do i = 0, npix1-1			# zero indexed
			Memr[oIm+i] = aimag (Memx[xline+i])
		    in = in + 1
		    if (in > ny)
			in = in - ny
		}
	    }
	    call sfree (sp)

	} else {				# just copy without centering

	    # Convert from complex, and write the data to output.
	    if (FT_REAL(fto) == YES && FT_IMAG(fto) == YES) {
		do j = 1, npix2 {
		    oRe = impl2r (oimRe, j)
		    oIm = impl2r (oimIm, j)
		    call aupxr (ximg[1,j], Memr[oRe], Memr[oIm], npix1)
		}
	    } else if (FT_REAL(fto) == YES) {
		do j = 1, npix2 {
		    oRe = impl2r (oimRe, j)
		    call achtxr (ximg[1,j], Memr[oRe], npix1)
		}
	    } else if (FT_IMAG(fto) == YES) {
		do j = 1, npix2 {
		    oIm = impl2r (oimIm, j)
		    do i = 1, npix1
			Memr[oIm+i-1] = aimag (ximg[i,j])
		}
	    }
	}
end

# ft_shift_x -- shift a complex 1-D array
# This routine copies an array from input to output, shifting by an
# integral number of pixels.  The actual arguments corresponding to
# the input and output arrays must be distinct.

procedure ft_shift_x (in, shift, out, npts)

complex in[npts]	# i: input array
int	shift		# i: amount of shift
complex out[npts]	# o: output array
int	npts		# i: size of arrays
#--
int	sh		# abs (shift)
int	n		# npts - abs(shift)

begin
	if (shift >= 0) {

	    n = npts - shift

	    call amovx (in, out[shift+1], n)
	    call amovx (in[n+1], out, shift)

	} else if (shift < 0) {

	    sh = abs (shift)
	    n = npts - sh

	    call amovx (in, out[n+1], sh)
	    call amovx (in[sh+1], out, n)
	}
end

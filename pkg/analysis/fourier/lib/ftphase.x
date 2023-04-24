include <imhdr.h>
include <mach.h>		# for EPSILOND
include <math.h>		# for TWOPI
include "../fourier.h"

# ft_phase_init -- initialize for phase shift
# The shift between the reference pixel and the first pixel is gotten.
# This is the shift for the forward transform; for the inverse transform
# change the sign of the shift for each axis.
#
# Phil Hodge, 16-Sep-1988  Subroutine created.
# Phil Hodge, 19-Feb-1990  Use mwcs instead of xt_wcs.

procedure ft_phase_init (ft, shift, not_zero)

pointer ft		# i: pointer to internal Fourier structure
double	shift[ARB]	# o: 1. - crpix
bool	not_zero[ARB]	# o: phase shift non-trivial?
#--
pointer im		# pointer to image header structure
int	k		# loop index for axis number

begin
	im = FT_IMAGE(ft)

	do k = 1, IM_NDIM(im) {

	    shift[k] = 1. - FT_CRPIX(ft,k)

	    if (abs (shift[k]) > EPSILOND)
		not_zero[k] = true		# shift is non-zero
	    else
		not_zero[k] = false
	}
end


# ft_phase_init_x -- initialize for cross correlation

procedure ft_phase_init_x (ft, ref_crval, ref_crpix,
		shift, not_zero)

pointer ft		# i: pointer to ft struct for first image
double	ref_crval[ARB]	# i: world coords at ref pixel in second image
double	ref_crpix[ARB]	# i: reference pixel in second image
double	shift[ARB]	# o: shift (e.g. 1-crpix) in pixels
bool	not_zero[ARB]	# o: phase shift non-trivial?
#--
pointer im		# pointer to imhdr struct for first image
pointer mw, ct		# pointers to world coordinate system struct
double	pix[2]		# pixel coordinates of ref_crval in first image
int	naxis		# dimension of image
int	k		# loop index for axis number
pointer mw_openim(), mw_sctran()
double	mw_c1trand()

begin
	im = FT_IMAGE(ft)

	naxis = IM_NDIM(im)

	# ref_crval is the set of coordinates at the reference pixel in
	# the second image; we want the pixel coordinates in the first
	# image of the point that has world coordinates ref_crval.
	mw = mw_openim (im)

	# Convert to pixel coordinates.
	if (naxis == 1) {
	    ct = mw_sctran (mw, "world", "logical", 1B)
	    pix[1] = mw_c1trand (ct, ref_crval[1])
	} else {
	    ct = mw_sctran (mw, "world", "logical", 3B)
	    call mw_c2trand (ct, ref_crval[1], ref_crval[2], pix[1], pix[2])
	}

	call mw_ctfree (ct)
	call mw_close (mw)

	do k = 1, naxis {

	    shift[k] = ref_crpix[k] - pix[k]

	    if (abs (shift[k]) > EPSILOND)
		not_zero[k] = true		# shift is non-zero
	    else
		not_zero[k] = false
	}
end


# ft_phase_tab -- fill cosine and sine tables
# This routine computes the cosines and sines for shifting an image
# by adjusting the phase of its Fourier transform and stores those
# values in arrays in the work struct.  This is called for each axis.

procedure ft_phase_tab (work, shift, not_zero)

pointer work		# i: structure for work space
double	shift		# i: shift (e.g. 1-crpix) for current axis
bool	not_zero	# i: phase shift for current axis non-trivial?
#--
pointer ct, st		# pointers to cosine, sine tables (in work struct)
double	xdimen		# =npts
double	phase_increment # change in phase for one pixel
double	phase		# phase at a particular pixel
int	npts		# length of current axis
int	hnpts		# npts / 2 + 1
int	k		# loop index

begin
	npts = NPTS(work)
	xdimen = double(npts)
	hnpts = npts / 2 + 1

	# Fill cosine and sine tables.
	phase_increment = TWOPI * shift / xdimen
	ct = COSTAB(work)
	st = SINTAB(work)
	Memr[ct] = 1.			# cos & sin for i = 1
	Memr[st] = 0.
	if (not_zero) {
	    do k = 1, hnpts-1 {		# do i = 2, hnpts	<-- NOTE!
		phase = phase_increment * k
		Memr[ct+k] = cos (phase)
		Memr[st+k] = sin (phase)
	    }
	} else {			# shift is zero
	    do k = 1, hnpts-1 {
		Memr[ct+k] = 1.
		Memr[st+k] = 0.
	    }
	}
end


# ft_phase_s -- shift by specified amount
#
# Apply phase shift in order to put coordinate zero at the first pixel.
#
# The shift theorem says that the transforms of f(x) = F(s) and
# f(x-a) = F(u) are related by :
#
#     F(u)  = exp (- 2 pi i a s) * F(s)
#           = [cos (2 pi a s) - i sin (2 pi a s)] * F(s)
#
# and since we are carting the real and imaginary parts around separately,
# we can compute the shifted transform :
#
# Re {F(u)} =   Re {F(s)} * cos (2 pi a s) + Im {F(s)} * sin (2 pi a s)
# Im {F(u)} = - Re {F(s)} * sin (2 pi a s) + Im {F(s)} * cos (2 pi a s)

procedure ft_phase_s (npts, realpt, imagpt, costab, sintab)

int	npts		# i: size of arrays
real	realpt[npts]	# io: real part
real	imagpt[npts]	# io: imaginary part
real	costab[ARB]	# i: cosine table (size = half of npts)
real	sintab[ARB]	# i: sine table (size = half of npts)
#--
real	w		# scratch
int	k		# loop index
int	j		# index of mirror image of k

begin
	# First half.
	do k = 2, npts/2 + 1 {
	    w         =   realpt[k] * costab[k] + imagpt[k] * sintab[k]
	    imagpt[k] = - realpt[k] * sintab[k] + imagpt[k] * costab[k]
	    realpt[k] = w
	}

	# Second half; use complex conjugate of (cos+i*sin).
	j = 2
	do k = npts, npts/2 + 2, -1 {
	    w         = realpt[k] * costab[j] - imagpt[k] * sintab[j]
	    imagpt[k] = realpt[k] * sintab[j] + imagpt[k] * costab[j]
	    realpt[k] = w
	    j = j + 1
	}
end

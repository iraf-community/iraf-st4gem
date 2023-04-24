include	<imhdr.h>
include "../fourier.h"

# These routines change the coordinate parameters (in the ft struct)
# for a forward or inverse Fourier transform.

# Phil Hodge, 23-Sep-1991	Modify to use MWCS.

# ft_wcs_fwd -- forward transform world coordinate system

# The reference pixel in the Fourier domain is set to one, and the
# coordinate value at the reference pixel is set to zero.  The coordinate
# type is transformed if appropriate (e.g. LAMBDA --> WAVENUMB).
# For a one-dimensional image, the pixel spacing in the Fourier domain
# is 1 / (N1 * CD1_1), where N1 is the length of the image.  For a two-
# dimensional image, let D1 and D2 be the pixel spacing in physical units
# for the two axes.  (I'm not sure what to do with different units in
# the two axes if the rotation is non-zero.)  Then:
# D1 = sqrt (CD1_1 ** 2 + CD2_1 ** 2)
# D2 = sqrt (CD1_2 ** 2 + CD2_2 ** 2)
# The new CD matrix is the old CD matrix multiplied on the right by the matrix

#  | 1/(N1*D1**2)     0         |
#  |    0          1/(N2*D2**2) |

procedure ft_wcs_fwd (fti, fto, ctd)

pointer fti		# i: FFT pointer; input image
pointer fto		# i: FFT pointer; output image
pointer ctd[2]		# i: pointers to input & output CTYPE dictionaries
#--
pointer im		# ptr to imhdr struct for input image
double	d[2]		# im_len * cdelt**2 for each axis
int	ndim, i
bool	fwd		# forward transform?  (true)

begin
	fwd = true
	im = FT_IMAGE(fti)
	ndim = IM_NDIM(im)

	do i = 1, ndim {
	    # crpix is changed later if center = yes.
	    FT_CRVAL(fto,i) = 0.d0
	    FT_CRPIX(fto,i) = 1.d0
	    call ft_change_ctype (fwd, Memc[ctd[1]], Memc[ctd[2]],
			FT_CTYPE(fti,i), FT_CTYPE(fto,i))
	}

	# Update the CD matrix.
	if (ndim == 1) {
	    FT_CD(fto,1,1) = 1.d0 / (FT_CD(fti,1,1) * IM_LEN(im,1))

	} else {			# 2-D image

	    d[1] = IM_LEN(im,1) * (FT_CD(fti,1,1) ** 2 + FT_CD(fti,2,1) ** 2)
	    d[2] = IM_LEN(im,2) * (FT_CD(fti,1,2) ** 2 + FT_CD(fti,2,2) ** 2)

	    FT_CD(fto,1,1) = FT_CD(fti,1,1) / d[1]
	    FT_CD(fto,1,2) = FT_CD(fti,1,2) / d[2]
	    FT_CD(fto,2,1) = FT_CD(fti,2,1) / d[1]
	    FT_CD(fto,2,2) = FT_CD(fti,2,2) / d[2]
	}
end

# ft_wcs_inv -- inverse transform world coordinate system

procedure ft_wcs_inv (fti, fto, ctd, coord_shift)

pointer fti		# i: FFT pointer; input image
pointer fto		# i: FFT pointer; output image
pointer ctd[2]		# i: pointers to input & output CTYPE dictionaries
bool	coord_shift	# io: set to false if ocrval, ocrpix not found
#--
pointer im		# ptr to imhdr struct for input image
double	d[2]		# im_len * cdelt**2 for each axis
int	ndim, i
bool	fwd		# forward transform? (false)
bool	foundit		# true if ocrval, ocrpix found in input file

begin
	fwd = false
	im = FT_IMAGE(fti)
	ndim = IM_NDIM(im)

	# Get original values of crpix & crval from keywords ocrpix & ocrval
	# in input image and assign to elements in fto structure.
	call load_old_info (FT_IMAGE(fti), fto, foundit)
	if ( coord_shift && !foundit ) {
	    call eprintf (
	"warning:  ocrpix & ocrval not found, so coord_shift reset to false\n")
	    coord_shift = false
	}

	# Update ctype.
	do i = 1, ndim
	    call ft_change_ctype (fwd, Memc[ctd[1]], Memc[ctd[2]],
			FT_CTYPE(fti,i), FT_CTYPE(fto,i))

	# Update the CD matrix.
	if (ndim == 1) {
	    FT_CD(fto,1,1) = 1.d0 / (FT_CD(fti,1,1) * IM_LEN(im,1))

	} else {			# 2-D image

	    d[1] = IM_LEN(im,1) * (FT_CD(fti,1,1) ** 2 + FT_CD(fti,2,1) ** 2)
	    d[2] = IM_LEN(im,2) * (FT_CD(fti,1,2) ** 2 + FT_CD(fti,2,2) ** 2)

	    FT_CD(fto,1,1) = FT_CD(fti,1,1) / d[1]
	    FT_CD(fto,1,2) = FT_CD(fti,1,2) / d[2]
	    FT_CD(fto,2,1) = FT_CD(fti,2,1) / d[1]
	    FT_CD(fto,2,2) = FT_CD(fti,2,2) / d[2]
	}
end

# ft_wcs_n -- assign world coordinate system
# Several of the coordinate parameters are copied from input to output,
# while crpix and crval are set to 1 and 0 respectively.  This is appropriate
# for powerspec or convolution.

procedure ft_wcs_n (fti, fto)

pointer fti		# i: FFT pointer; input image
pointer fto		# i: FFT pointer; output image
#--
pointer im		# ptr to imhdr struct for input image
int	ndim, i, j

begin
	im = FT_IMAGE(fti)
	ndim = IM_NDIM(im)

	do i = 1, ndim {
	    FT_CRVAL(fto,i) = 0.d0
	    FT_CRPIX(fto,i) = 1.d0
	    call strcpy (FT_CTYPE(fti,i), FT_CTYPE(fto,i), SZ_CTYPE)
	    do j = 1, ndim
		FT_CD(fto,i,j) = FT_CD(fti,i,j)
	}
end

# ft_wcs_0 -- assign default world coordinate system
# The coordinate parameters are set to zero or one.
# This is appropriate for cross correlation or autocorrelation.

procedure ft_wcs_0 (ft)

pointer ft		# i: FFT pointer
#--
pointer im		# ptr to imhdr struct for image
int	ndim, i, j

begin
	im = FT_IMAGE(ft)
	ndim = IM_NDIM(im)

	do i = 1, ndim {
	    FT_CRVAL(ft,i) = 0.d0
	    FT_CRPIX(ft,i) = 1.d0
	    call strcpy ("PIXEL", FT_CTYPE(ft,i), SZ_CTYPE)
	    do j = 1, ndim {
		if (i == j)
		    FT_CD(ft,i,j) = 1.d0
		else
		    FT_CD(ft,i,j) = 0.d0
	    }
	}
end

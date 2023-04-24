# This file contains routines for getting coordinate info from an image
# into the ft struct or putting info from ft to an image.  The following
# routines are included:
#
#  load_ctstruct - from image into ft struct
#  save_ctstruct - from ft struct into image header
#  load_old_info - copy ocrpix & ocrval from image into ft struct
#  save_old_info - copy crpix & crval to ocrpix & ocrval in image header
#
# C. Biemesderfer		Subroutines created.
# Phil hodge, 14-Nov-1990	Check for cdelt=0 for 1-D, and print
#				warning instead of calling error.
# Phil Hodge, 23-Sep-1991	Modify to use MWCS.
# Phil Hodge, 22-Nov-1991	In load_ctstruct, use imaccf for CTYPEi.
# Phil Hodge, 19-May-1992	Swap cd[1,2] and cd[2,1] when copying between
#				mwcs variables and FT_CD because of different
#				convention for storing matrices.  Also, in
#				save_ctstruct, use imadd instead of mw_saveim
#				to save (or add) the coordinate parameters.
# Phil Hodge, 24-Sep-1992	In load_ctstruct use wcsdim and call
#				ft_ct_extract; simplify save_ctstruct.
# Phil Hodge, 16-Oct-1992	In ft_ct_extract assign default axis numbers
#				if values are not found; this is a temporary
#				fix for running on a vax.
# Phil Hodge, 27-Jun-1994	Call stx_getcoord in load_ctstruct; don't use
#				mw_saveim with null mwcs in save_ctstruct.

include	<imhdr.h>
include	<math.h>
include "../fourier.h"

define	SZ_PNAME	8

# load_ctstruct -- load coordinate information
# Get the size and coordinate information from the image, and load
# that info into the ft structure.

procedure load_ctstruct (ft, im)

pointer ft			# i: pointer to ft structure
pointer im			# i: pointer to image
#--
double	crpix[IM_MAXDIM]	# reference pixel
double	crval[IM_MAXDIM]	# world coordinates at reference pixel
double	cd[IM_MAXDIM,IM_MAXDIM]	# CD matrix
char	ctype[SZ_CTYPE,IM_MAXDIM]	# coordinate type
int	ndim			# dimension of image (or image section)
int	i, j
errchk	stx_getcoord

begin
	ndim = IM_NDIM(im)
	FT_NAXIS(ft) = ndim

	# Get the coordinate parameters for this image.  stx_getcoord
	# adjusts these parameters, if necessary, to account for any
	# image section that was used with the image name.
	call stx_getcoord (im, crpix, crval, cd, IM_MAXDIM, ctype, SZ_CTYPE)

	# Copy info to the ft structure.
	do i = 1, ndim {
	    FT_CRPIX(ft,i) = crpix[i]
	    FT_CRVAL(ft,i) = crval[i]
	    call strcpy (ctype[1,i], FT_CTYPE(ft,i), SZ_CTYPE)
	    do j = 1, ndim
		FT_CD(ft,i,j) = cd[i,j]
	}
end

# SAVE_CTSTRUCT -- Copy coordinate information from ft structure into image.

procedure save_ctstruct (ft, im)

pointer ft	# i: pointer to ft structure
pointer im	# i: pointer to image
#--
char	pname[SZ_PNAME]
int	ndim
int	wcsdim
int	i
int	imaccf(), imgeti()

begin
	ndim = IM_NDIM(im)

	# Reset some of the MWCS info.  (Leave LTM, LTV, WAT keywords.)

	# Set WCS dimension to the correct value.
	if (imaccf (im, "wcsdim") == YES) {
	    wcsdim = imgeti (im, "wcsdim")
	    call imputi (im, "wcsdim", ndim)
	} else {
	    wcsdim = ndim
	}

	# Remove any WAXMAP keywords.
	do i = 1, 99 {
	    call sprintf (pname, SZ_PNAME, "WAXMAP%02d")
		call pargi (i)
	    if (imaccf (im, pname) == YES)
		call imdelf (im, pname)
	    else
		break
	}

	# Save coordinate parameters in the header.
	if (ndim == 1) {
	    call imaddd (im, "crpix1", FT_CRPIX(ft,1))
	    call imaddd (im, "crval1", FT_CRVAL(ft,1))
	    call imaddd (im, "cd1_1", FT_CD(ft,1,1))
	    call imastr (im, "ctype1", FT_CTYPE(ft,1))
	} else if (ndim == 2) {
	    call imaddd (im, "crpix1", FT_CRPIX(ft,1))
	    call imaddd (im, "crpix2", FT_CRPIX(ft,2))
	    call imaddd (im, "crval1", FT_CRVAL(ft,1))
	    call imaddd (im, "crval2", FT_CRVAL(ft,2))
	    call imaddd (im, "cd1_1", FT_CD(ft,1,1))
	    call imaddd (im, "cd1_2", FT_CD(ft,1,2))
	    call imaddd (im, "cd2_1", FT_CD(ft,2,1))
	    call imaddd (im, "cd2_2", FT_CD(ft,2,2))
	    call imastr (im, "ctype1", FT_CTYPE(ft,1))
	    call imastr (im, "ctype2", FT_CTYPE(ft,2))
	}
end


# load_old_info -- load old coordinate info
# This routine gets the original values of crval & crpix from the input
# image by the keywords ocrval & ocrpix (with the axis number appended,
# of course) and loads the values into the FT structure.

procedure load_old_info (im, fto, foundit)

pointer im		# i: pointer to output image
pointer fto		# i: ft structure for output image
bool	foundit		# o: true if both ocrval & ocrpix were found
#--
char	pname[SZ_PNAME]
int	nax, iax
double	imgetd()
int	imaccf()

begin
	foundit = true			# initial value

	nax = IM_NDIM(im)

	do iax = 1, nax {

	    # CRVAL for each axis
	    call sprintf (pname, SZ_PNAME, "ocrval%d")
		call pargi (iax)
	    if (imaccf (im, pname) == YES) {
		FT_CRVAL(fto,iax) = imgetd (im, pname)
	    } else {
		FT_CRVAL(fto,iax) = 0.d0
		foundit = false
	    }

	    # CRPIX for each axis
	    call sprintf (pname, SZ_PNAME, "ocrpix%d")
		call pargi (iax)
	    if (imaccf (im, pname) == YES) {
		FT_CRPIX(fto,iax) = imgetd (im, pname)
	    } else {
		FT_CRPIX(fto,iax) = 1.d0
		foundit = false
	    }
	}
end

# save_old_info -- save old coordinate info
# This routine copies the values of crval & crpix from the input image
# to the output image but with keywords ocrval & ocrpix (with the axis
# number appended, of course).

procedure save_old_info (fti, im)

pointer fti	# i: ft structure for input image
pointer im	# i: pointer to output image
#--
char	pname[SZ_PNAME]
int	nax, iax

begin
	nax = IM_NDIM(im)

	do iax = 1, nax {

	    # CRVAL for each axis
	    call sprintf (pname, SZ_PNAME, "ocrval%d")
		call pargi (iax)
	    call imaddd (im, pname, FT_CRVAL(fti,iax))

	    # CRPIX for each axis
	    call sprintf (pname, SZ_PNAME, "ocrpix%d")
		call pargi (iax)
	    call imaddd (im, pname, FT_CRPIX(fti,iax))
	}
end

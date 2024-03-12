###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	procedure getimage (image)
#		char	image[SZ_FNAME]
#
#		procedure getval(line, rval)
#		char	line[ARB]
#		real	rval
#
#  Description:	GETIMAGE -- Opens IRAF images for gettemps and getspec
#		and returns the image file descriptors and extracted
#		header information.
#		Adapted from the onedspec package in IRAF.
#
#		GETVAL -- returns a real value from an ASCII header line
#		Adapted from the onedspec package in IRAF.
#
#  Arguments:	
#
#  Returns:	
#
#  Notes:
#
#  History:	June	1987	Gerard Kriss
#		Oct  1  1991	gak	Adapted for specfit.  Handles 1D and
#					2D (HUT ballistic output) IRAF images
#		June 16	1995	grimes	increased header size
#
###########################################################################

include	<imhdr.h>
include "specfit.h"

define	COL_VALUE	11	#Start column for FITS data
define	LEN_KEY		8	#Length of FITS keyword
define	USER_AREA	Memc[($1+IMU-1)*SZ_STRUCT + 1]
define	LEN_USER_AREA	8640

# GETIMAGE -- Read new image pixels

procedure getimage (image)
char	image[SZ_FNAME]

char	line[SZ_LINE]
int	fd, nline, i
real	w0, wpc

int	stropen(), getline(), strncmp()
real	abs(), sqrt()
pointer	immap(), imgl1r(), imgl2r()
pointer im, fptr, eptr

include "specfit.com"

begin
	# Map the image
	iferr (im = immap (image, READ_ONLY, LEN_USER_AREA)) {
	    call eprintf ("Cannot open image: %s\n")
	    call pargstr (image)
	    im = ERR
	    return
	}

	# Get header info
	npts = IM_LEN (im, 1)
	fd = stropen (USER_AREA(im), (LEN_USER_AREA - 1)*SZ_STRUCT, READ_ONLY)

	# Allocate storage
	call malloc(spectrum, npts, TY_REAL)
	call malloc(lambdas, npts, TY_REAL)
	call malloc(errors, npts, TY_REAL)

	#Initialize defaults
	w0 = 0.
	wpc = 1.
	call strcpy (IM_TITLE(im), specname, SZ_FNAME)

	# Get keywords
	while (getline (fd, line) != EOF) {
	    if (strncmp (line, "W0      ", LEN_KEY) == 0) {
		call getval (line, w0)
	    } else if (strncmp (line, "WPC     ", LEN_KEY) == 0) {
		call getval (line, wpc)
	    } else if (strncmp (line, "CRVAL1  ", LEN_KEY) == 0) {
		call getval (line, w0)
	    } else if (strncmp (line, "CDELT1  ", LEN_KEY) == 0) {
		call getval (line, wpc)
	    }
	}

	# Get fluxes, wavelengths, and errors
	if (IM_NDIM(im) == 2) {
	    nline = 1
	    fptr = imgl2r (im, nline)
	    nline = 2
	    eptr = imgl2r (im, nline)
	} else if (IM_NDIM(im) == 1) {
	    nline = 1
	    fptr = imgl1r (im, 1)
	} else {
		call eprintf("Image %s has too many dimensions (%d)\n")
			call pargstr(image)
			call pargi(IM_NDIM(im))
		return
	}

	if ( nline == 1 ) {
		for ( i = 0; i < npts; i = i + 1 ) {
			Memr[lambdas+i] = w0 + i * wpc
			Memr[spectrum+i] = Memr[fptr+i]
			Memr[errors+i] = sqrt(v0 + v1 * abs(Memr[fptr+i]))
		}
	} else {
		for ( i = 0; i < npts; i = i + 1 ) {
			Memr[lambdas+i] = w0 + i * wpc
			Memr[spectrum+i] = Memr[fptr+i]
			Memr[errors+i] = Memr[eptr+i]
		}
	}

	# Unmap the image
	call imunmap(im)

end

procedure getval(line, rval)

char	line[ARB]
real	rval

int	ip
double	dval

int	ctod()

begin
	ip = COL_VALUE
	if ( ctod(line, ip, dval) != 0 )
		rval = dval
	else
		rval = INDEFR
end

include <imhdr.h>
include	"simtwo.h"
include "function.h"

#* HISTORY *
#* B.Simon	04-May-95	Original
#* B.Simon	03-Oct-95	Test for zero pixel size

# FNIMAGE.X -- Functions for using extended objects read from images

# DELIMAGE -- Delete variables used by image object

procedure delimage ()

#--
include	"fnimage.com"

begin
	call imunmap (im)
	call mfree (grid, TY_REAL)
end

# NEWIMAGE -- Initialize variables used by image object

procedure newimage (func, apscale, dynrange, nsub)

pointer	func		# i: function descriptor
double	apscale		# i: aperture scale
real	dynrange	# i: dynamic range of object
int	nsub		# i: number of subpixels
#--
include	"fnimage.com"

int	ixlo, ixhi, iylo, iyhi
pointer	sp, image, object
real	xpcen, ypcen, range, imscale, gscale

string	radform  RADIUS_UNITS
string	badargs  "Argument type/number mismatch"
string	badsize  "Image has invalid scale; detector scale assumed"

int	realargs()
pointer	immap(), imgs2r()

begin
	# Allocate string for reconstructing function name

	call smark (sp)
	call salloc (image, SZ_FNAME, TY_CHAR)
	call salloc (object, SZ_LINE, TY_CHAR)

	# Extract and validate arguments from function descriptor

	if (FUN_TYPE(func,2) != TY_CHAR) {
	    call namefunc (func, Memc[object], SZ_LINE)
	    call printerr_str (badargs, Memc[object])
	}

	if (realargs (func, 3, radius, 1) == NO) {
	    call printerr_str (badargs, Memc[object])
	    call namefunc (func, Memc[object], SZ_LINE)
	}
	call angtodegr (radform, radius)

	# Open image

	im = immap (FUN_STR(func,2), READ_ONLY, 0)

	# Compute interpolation grid

	call pixregion (im, dynrange, xpcen, ypcen, npix)
	call pixsize (im, imscale)

	if (imscale <= 0.0) {
	    call strcpy (IM_HDRFILE(im), Memc[image], SZ_FNAME)
	    call synphotwarn (badsize, Memc[image])
	    imscale = apscale
	}

	range = imscale * npix

	npix = nsub * (int (radius / apscale) + 1)

	call malloc (grid, npix+1, TY_REAL)

	imscale = imscale * radius / range
	gscale = apscale / (nsub * imscale)

	call setgrid (gscale, Memr[grid], npix+1)

	# Read buffer containing image

	ixlo = max (int (xpcen + Memr[grid]), 1)
	ixhi = min (int (xpcen + Memr[grid+npix] + 1.0), IM_LEN(im,1))
	iylo = max (int (ypcen + Memr[grid]), 1)
	iyhi = min (int (ypcen + Memr[grid+npix] + 1.0), IM_LEN(im,2))
	nx = (ixhi - ixlo) + 1
	ny = (iyhi - iylo) + 1

	buf = imgs2r (im, ixlo, ixhi, iylo, iyhi)

	# Initialize image center

	xcen = xpcen - (ixlo - 1)
	ycen = ypcen - (iylo - 1)

	call sfree (sp)
end

# PNTIMAGE -- Return the value of a image at a point

real procedure pntimage (x, y)

real	x		# i: x coordinate
real	y		# i: y coordinate
#--
include "fnimage.com"

int	ix, iy, p1, p2, p3, p4
real	ax, bx, ay, by, z

begin
	# Locate grid points which bound the interpolated point

	ix = min (nx-1, max (1, int(x)))
	iy = min (ny-1, max (1, int(y)))

	p1 = ny * (iy - 1) + (ix - 1)
	p3 = p1 + ny
	p2 = p1 + 1
	p4 = p3 + 1

	# Compute interpolated value

	bx = x - ix
	ax = 1.0 - bx

	by = y - iy
	ay = 1.0 - by

	z = ax * ay * Memr[buf+p1] + bx * ay * Memr[buf+p2] +
	    ax * by * Memr[buf+p3] + bx * by * Memr[buf+p4]

	return (max (z, 0.0))
end

# SIZIMAGE -- Return size of image object

procedure sizimage (nox, noy)

int	nox		# o: number of pixels on x dimension of object
int	noy		# o: number of pixels on y dimension of object
#--
include "fnimage.com"

begin
	nox = npix
	noy = npix
end

# VALIMAGE -- Return value of pixel in image object

procedure valimage (ix, iy, sum, val)

int	ix		# i: x pixel index
int	iy		# i: y pixel index
real	sum		# u: partial sum of pixel value
real	val		# o: pixel value
#--
include "fnimage.com"

real	xlo, xhi, ylo, yhi

real    pntimage() # Return the value of a image at a point
extern	pntimage

begin
	# Compute endpoints of pixel

	xlo = Memr[grid+ix-1] + xcen
	xhi = Memr[grid+ix] + xcen

	ylo = Memr[grid+iy-1] + ycen
	yhi = Memr[grid+iy] + ycen

	# Initialize computation

	if (ix == 1)
	    sum = INDEFR

	# Integrate function's value over pixel

	call simpson2 (xlo, xhi, ylo, yhi, pntimage, sum, val)
end

include	"simtwo.h"

#* HISTORY *
#* B.Simon	04-May-95	Original

# FNGAUSS.X -- Functions for using extended objects with gaussian profiles

# DELGAUSS -- Delete variables used by gaussian object

procedure delgauss ()

#--
include	"fngauss.com"

begin
	call mfree (xgrid, TY_REAL)
	call mfree (ygrid, TY_REAL)
end

# NEWGAUSS -- Initialize variables used by gaussian object

procedure newgauss (shape, apscale, dynrange, nsub)

pointer	shape		# i: shape descriptor
double	apscale		# i: aperture scale
real	dynrange	# i: dynamic range of object
int	nsub		# i: number of subpixels
#--
include	"fngauss.com"

int	npix
pointer	sp, object
real	gscale, arglist[3]

string	radform  RADIUS_UNITS
string	badargs  "Argument type/number mismatch"

int	realargs()

begin
	# Allocate string for reconstructing shape name

	call smark (sp)
	call salloc (object, SZ_LINE, TY_CHAR)

	# Extract and validate arguments from shape descriptor

	if (realargs (shape, 2, arglist, 3) == NO) {
	    call namefunc (shape, Memc[object], SZ_LINE)
	    call printerr_str (badargs, Memc[object])
	}

	# Set common block variables

	radius = arglist[1]
	call angtodegr (radform, radius)
	ar = arglist[2]
	pa = arglist[3]
	scale = sqrt (log (2.0)) * (apscale / radius)

	# Compute object grid points

	npix = nsub * (int (2.0 * sqrt (log (dynrange)) / scale) + 1)
	nx = npix
	ny = npix

	call malloc (xgrid, nx+1, TY_REAL)
	call malloc (ygrid, ny+1, TY_REAL)

	gscale =  1.0 / nsub
	call setgrid (gscale, Memr[xgrid], nx+1)
	call setgrid (gscale, Memr[ygrid], ny+1)

	call sfree (sp)
end

# PNTGAUSS -- Return the value of a gaussian at a point

real procedure pntgauss (x, y)

real	x		# i: x coordinate
real	y		# i: y coordinate
#--
include "fngauss.com"

real	r, z

begin
	call radcalc (x, y, r)
	z = exp (- ((r * scale) ** 2))
	return (z)

end

# SIZGAUSS -- Return size of gaussian object

procedure sizgauss (nox, noy)

int	nox		# o: number of pixels on x dimension of object
int	noy		# o: number of pixels on y dimension of object
#--
include "fngauss.com"

begin
	nox = nx
	noy = ny
end

# VALGAUSS -- Return value of pixel in gaussian object

procedure valgauss (ix, iy, sum, val)

int	ix		# i: x pixel index
int	iy		# i: y pixel index
real	sum		# u: partial sum of pixel value
real	val		# o: pixel value
#--
include "fngauss.com"

real	xlo, xhi, ylo, yhi

real    pntgauss() # Return the value of a gaussian at a point
extern	pntgauss

begin
	# Compute endpoints of pixel

	xlo = Memr[xgrid+ix-1]
	xhi = Memr[xgrid+ix]

	ylo = Memr[ygrid+iy-1]
	yhi = Memr[ygrid+iy]

	# Initialize computation

	if (ix == 1) {
	    sum = INDEFR
	    if (iy == 1)
		call initrad (ar, pa)
	}

	# Integrate shape's value over pixel

	call simpson2 (xlo, xhi, ylo, yhi, pntgauss, sum, val)
end

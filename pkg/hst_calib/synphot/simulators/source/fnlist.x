include	<tbset.h>
include	"simtwo.h"
include "function.h"

define	NPROF		1001

#* HISTORY *
#* B.Simon	04-May-95	Original

# FNLIST.X -- Functions for using extended objects with a list profile

# DELLIST -- Delete variables used by list profile object

procedure dellist ()

#--
include	"fnlist.com"

begin
	call mfree (prof, TY_REAL)
	call mfree (xgrid, TY_REAL)
	call mfree (ygrid, TY_REAL)
end

# NEWLIST -- Initialize variables used by list profile object

procedure newlist (func, apscale, dynrange, nsub)

pointer	func		# i: function descriptor
double	apscale		# i: aperture scale
real	dynrange	# i: dynamic range of object
int	nsub		# i: number of subpixels
#--
include	"fnlist.com"

int	nrow, irow, npix
pointer	tp, cp, sp, object, table, xold, yold, y2, xnew, nulflg
real	gscale, arglist[3]

string	radform  RADIUS_UNITS
string	badargs  "Argument type/number mismatch"
string	listcol  "profile"
string	nullval  "Indef found in profile table"

int	realargs(), tbpsta()
pointer	opnsyntab()

begin
	# Allocate string for reconstructing function name

	call smark (sp)
	call salloc (object, SZ_LINE, TY_CHAR)
	call salloc (table, SZ_FNAME, TY_CHAR)

	# Extract and validate arguments from function descriptor

	if (FUN_TYPE(func,2) != TY_CHAR) {
	    call namefunc (func, Memc[object], SZ_LINE)
	    call printerr_str (badargs, Memc[object])
	}

	call strcpy (FUN_STR(func,2), Memc[table], SZ_FNAME)

	if (realargs (func, 3, arglist, 3) == NO) {
	    call namefunc (func, Memc[object], SZ_LINE)
	    call printerr_str (badargs, Memc[object])
	}

	# Set common block variables

	radius = arglist[1]
	call angtodegr (radform, radius)
	ar = arglist[2]
	pa = arglist[3]
	coef = radius / apscale
	scale = (NPROF - 1) * apscale / radius

	# Open profile table

	tp = opnsyntab (Memc[table])
	nrow = tbpsta (tp, TBL_NROWS)

	# Allocate memory for arrays

	call salloc (xold, nrow, TY_REAL)
	call salloc (yold, nrow, TY_REAL)
	call salloc (y2, nrow, TY_REAL)
	call salloc (nulflg, nrow, TY_BOOL)
	call salloc (xnew, NPROF+1, TY_REAL)
	call malloc (prof, NPROF+1, TY_REAL)

	# Read profiles from table

	call syncolptr (tp, listcol, 1, cp)
	call tbcgtr (tp, cp, Memr[yold], Memb[nulflg], 1, nrow)

	do irow = 1, nrow {
	    if (Memb[nulflg+irow-1])
		call printerr_str (nullval, Memc[table])
	}

	# Use spline interpolation to compute profile

	do irow = 0, nrow-1
	    Memr[xold+irow] = irow

	do irow = 0, NPROF-1
	    Memr[xnew+irow] = real (irow * (nrow - 1)) / (NPROF - 1)

	call synspline (nrow, Memr[xold], Memr[yold], Memr[y2])

	call evalspline (nrow, Memr[xold], Memr[yold], Memr[y2], 
			 NPROF, Memr[xnew], Memr[prof])

	Memr[prof+NPROF] = 0.0

	# Compute object grid points

	npix = nsub * (int (2.0 * radius / apscale) + 1)
	nx = npix
	ny = npix

	call malloc (xgrid, nx+1, TY_REAL)
	call malloc (ygrid, ny+1, TY_REAL)

	gscale =  1.0 / nsub
	call setgrid (gscale, Memr[xgrid], nx+1)
	call setgrid (gscale, Memr[ygrid], ny+1)

	call sfree (sp)
end

# PNTLIST -- Return the value of a list profile at a point

real procedure pntlist (x, y)

real	x		# i: x coordinate
real	y		# i: y coordinate
#--
include "fnlist.com"

int	idx
real	r, z, frac, rs

begin
	call radcalc (x, y, r)

	# Use linear interpolation to calculate value

	if (r > coef) {
	    z = 0.0
	} else {
	    rs = r * scale
	    idx = rs
	    frac = rs - idx
	    z = Memr[prof+idx] * (1.0 - frac) + Memr[prof+idx+1] * frac
	}

	return (z)

end

# SIZLIST -- Return size of list profile object

procedure sizlist (nox, noy)

int	nox		# o: number of pixels on x dimension of object
int	noy		# o: number of pixels on y dimension of object
#--
include "fnlist.com"

begin
	nox = nx
	noy = ny
end

# VALLIST -- Return value of pixel in list profile object

procedure vallist (ix, iy, sum, val)

int	ix		# i: x pixel index
int	iy		# i: y pixel index
real	sum		# u: partial sum of pixel value
real	val		# o: pixel value
#--
include "fnlist.com"

real	xlo, xhi, ylo, yhi

real    pntlist() # Return the value of a list profile at a point
extern	pntlist

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

	# Integrate function's value over pixel

	call simpson2 (xlo, xhi, ylo, yhi, pntlist, sum, val)
end

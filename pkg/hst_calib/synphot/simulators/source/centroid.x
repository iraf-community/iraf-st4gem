include <imhdr.h>

#* HISTORY *
#* B.Simon	28-Sep-95	Original

# CENTROID -- Compute star center using MPC algorithm.
#
# Code was shamelessly stolen from the imcntr task in the iraf 
# centroid package

procedure centroid2 (im, xstart, ystart, boxsize, xcntr, ycntr)

pointer	im		# i: image descriptor
real	xstart		# i: initial guess at x center
real	ystart		# i: initial guess at y center
int	boxsize		# i: box size used in centroiding
real	xcntr		# o: computed x center
real	ycntr		# o: computed y center
#--
int	x1, x2, y1, y2, half_box
int	ncols, nrows, nx, ny, try
pointer	bufptr, sp, x_vect, y_vect
real	xinit, yinit

int	imgs2r()

begin
	half_box = (boxsize - 1) / 2
	xinit = xstart
	yinit = ystart

	# Mark region to extract - use box size
	ncols = IM_LEN (im, 1)
	nrows = IM_LEN (im, 2)
	try = 0

	repeat {
	    x1 = amax1 (xinit - half_box, 1.0) +0.5
	    x2 = amin1 (xinit + half_box, real(ncols)) +0.5
	    y1 = amax1 (yinit - half_box, 1.0) +0.5
	    y2 = amin1 (yinit + half_box, real(nrows)) +0.5
	    nx = x2 - x1 + 1
	    ny = y2 - y1 + 1

	    # Extract region around center
	    bufptr = imgs2r (im, x1, x2, y1, y2)

	    # Collapse to two 1-D arrays
	    call smark (sp)
	    call salloc (x_vect, nx, TY_REAL)
	    call salloc (y_vect, ny, TY_REAL)

	    call aclrr (Memr[x_vect], nx)
	    call aclrr (Memr[y_vect], ny)

	    # Sum all rows
	    call cen_rowsum (Memr[bufptr], Memr[x_vect], nx, ny)

	    # Sum all columns
	    call cen_colsum (Memr[bufptr], Memr[y_vect], nx, ny)

	    # Find centers
	    call cen_getcenter (Memr[x_vect], nx, xcntr)
	    call cen_getcenter (Memr[y_vect], ny, ycntr)
	    call sfree (sp)

	    # Check for INDEF centers.
	    if (IS_INDEFR(xcntr) || IS_INDEFR(ycntr)) {
		xcntr = xinit
		ycntr = yinit
		break
	    }

	    # Add in offsets
	    xcntr = xcntr + x1
	    ycntr = ycntr + y1

	    try = try + 1
	    if (try == 1) {
		if ((abs(xcntr-xinit) > 1.0) || (abs(ycntr-yinit) > 1.0)) {
		    xinit = xcntr
		    yinit = ycntr
		}
	    } else
		break
	}
end

procedure centroid1 (im, start, boxsize, center)

pointer	im		# i: image descriptor
real	start		# i: initial guess at center
int	boxsize		# i: box size used in centroiding
real	center		# o: computed center
#--
int	half_box, npix, x1, x2, nx, try
pointer	bufptr
real	xinit

int	imgs1r()

begin
	half_box = (boxsize - 1) / 2
	xinit = start

	# Mark region to extract - use box size
	npix = IM_LEN (im, 1)
	try = 0

	repeat {
	    x1 = amax1 (xinit - half_box, 1.0) +0.5
	    x2 = amin1 (xinit + half_box, real(npix)) +0.5
	    nx = x2 - x1 + 1

	    # Extract region around center
	    bufptr = imgs1r (im, x1, x2)

	    # Find centers
	    call cen_getcenter (Memr[bufptr], nx, center)

	    # Check for INDEF centers.
	    if (IS_INDEFR(center)) {
		center = xinit
		break
	    }

	    # Add in offsets
	    center = center + x1

	    try = try + 1
	    if (try == 1) {
		if ((abs(center - xinit) > 1.0)) {
		    xinit = center
		}
	    } else
		break
	}
end


# CEN_ROWSUM -- Sum all rows in a raster

procedure cen_rowsum (v, row, nx, ny)

int	nx, ny
real	v[nx,ny]
real	row[ARB]

int	i, j

begin
	do i = 1, ny
	    do j = 1, nx
		row[j] = row[j] + v[j,i]
end


# CEN_COLSUM -- Sum all columns in a raster.

procedure cen_colsum (v, col, nx, ny)

int	nx, ny
real	v[nx,ny]
real	col[ARB]

int	i, j

begin
	do i = 1, ny
	    do j = 1, nx
		col[j] = col[j] + v[i,j]
end


# CEN_GETCENTER -- Compute center of gravity of array.

procedure cen_getcenter (v, nv, vc)

real	v[ARB]
int	nv
real	vc

int	i
real	sum1, sum2, sigma, cont

begin
	# Assume continuum level is at endpoints
	# Compute first moment
	sum1 = 0.0
	sum2 = 0.0

	call aavgr (v, nv, cont, sigma)

	do i = 1, nv
	    if (v[i] > cont) {
	        sum1 = sum1 + (i-1) * (v[i] - cont)
	        sum2 = sum2 + (v[i] - cont)
	    }

	# Determine center
	if (sum2 == 0.0)
	    vc = INDEFR
	else
	    vc = sum1 / sum2
end

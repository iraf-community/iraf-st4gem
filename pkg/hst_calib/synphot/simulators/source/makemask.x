include	<math.h>
include "simtwo.h"
include "function.h"

define	MAXARGS		20
define	RECT		1
define	BARRED		2
define	PLANETARY	3
define	MULTI		4

# MAKEMASK -- Construct an aperture mask for the object

int procedure makemask (apshape, apscale, ox, oy, nsub, 
			nix, niy, ntx, nty, mask)

pointer	apshape		# i: aperture shape descriptor
double	apscale		# i: detector pixel size
real	ox		# i: x coord of object center
real	oy		# i: y coord of object center
int	nsub		# i: number of subpixels
int	nix		# i: first dimension of output buffer
int	niy		# i: second dimension of output buffer 
int	ntx		# i: x length of convolved object
int	nty		# i: y length of convolved object
real	mask[ntx,nty]	# o: aperture mask
#--
int	index, iarg, jarg, nargs, flag, status
pointer	sp, func
real	xlen, ylen, xstart, ystart
real	slit[2,5], arglist[MAXARGS]

string	apform	  APER_UNITS
string	aperlist  "rect,barred,planet,multi"
string	toomany   "Too many arguments for aperture shape function"
string	badtype   "Non-numeric argument in aperture shape function"

int	realargs(), word_match()
real	asumr()

begin
	# Check for no aperture case

	if (apshape == NULL) {
	    call amovkr (1.0, mask, ntx*nty)
	    return (YES)
	}

	call amovkr (0.0, mask, ntx*nty)

	# Allocate string for reconstructing aperture shape function

	call smark (sp)
	call salloc (func, SZ_LINE, TY_CHAR)

	# Convert aperture shape parameters into pixel units

	nargs = FUN_NPAR(apshape) - 1

	if (nargs > MAXARGS) {
	    call namefunc (apshape, Memc[func], SZ_LINE)
	    call printerr_str (toomany, Memc[func])
	}

	call strfix (FUN_STR(apshape,1))
	index = word_match (FUN_STR(apshape,1), aperlist)

	if (realargs (apshape, 2, arglist, nargs) == NO) {
	    call namefunc (apshape, Memc[func], SZ_LINE)
	    call printerr_str (badtype, Memc[func])
	}

	# Don't convert last argument for planetary slits,
	# as it is the slit tilt angle

	if (index == PLANETARY) {
	    jarg = nargs - 1
	} else {
	    jarg = nargs
	}

	do iarg = 1, jarg {
	    call angtodegr (apform, arglist[iarg])
	    arglist[iarg] = arglist[iarg] / apscale
	}

	# Create list of rectangles from apertures and mask each one
	
	switch (index) {
	case RECT:
	    # arguments are slit width (y dimension) and length (x dimension)
	    xlen = arglist[2]
	    ylen = arglist[1]
	    xstart = -0.5 * xlen
	    ystart = -0.5 * ylen
	    call makebox (nix, niy, xstart, ystart, xlen, ylen, slit)

	    call polymask (5, slit, ox, oy, nsub, ntx, nty, mask, flag)
	    status = flag

	case BARRED:
	    # arguments are slit width (y dimension) and slit length
	    # alternating with gap length (x dimension) from left 
	    # to right

	    ylen = arglist[1]
	    xstart = -0.5 * asumr (arglist[2], nargs-1)
	    ystart = -0.5 * ylen

	    status = NO
	    do iarg = 2, nargs, 2 {
		xlen = arglist[iarg]
		call makebox (nix, niy, xstart, ystart, xlen, ylen, slit)

		call polymask (5, slit, ox, oy, nsub, ntx, nty, mask, flag)
		if (flag == YES)
		    status = YES

		if (iarg < nargs)
		    xstart = xstart + xlen + arglist[iarg+1]
	    }

	case PLANETARY:
	    # arguments are slit width and length for each rectangular
	    # piece from left to right and rotation angle in degrees

	    xstart = 0.0
	    do iarg = 2, nargs-1, 2 
		xstart = xstart + arglist[iarg]
	    xstart = -0.5 * xstart

	    do iarg = 1, nargs-1, 2 {
		xlen = arglist[iarg+1]
		ylen = arglist[iarg]
		ystart = -0.5 * ylen

		call makebox (nix, niy, xstart, ystart, xlen, ylen, slit)
		call rotatbox (nix, niy, slit, arglist[nargs])

		call polymask (5, slit, ox, oy, nsub, ntx, nty, mask, flag)
		if (flag == YES)
		    status = YES

		xstart = xstart + xlen
	    }
	    
	case MULTI:
	    # arguments are slit widths, lengths, y offset to midpoint 
	    # of each slit and x distance between each slit

	    xlen = arglist[2]
	    ylen = arglist[1]
	    xstart = xlen

	    do iarg = 4, nargs, 2 
		xstart = xstart + arglist[iarg]
	    xstart = -0.5 * xstart

	    do iarg = 3, nargs, 2 {
		ystart = arglist[iarg] - 0.5 * ylen
		call makebox (nix, niy, xstart, ystart, xlen, ylen, slit)

		call polymask (5, slit, ox, oy, nsub, ntx, nty, mask, flag)
		if (flag == YES)
		    status = YES

		if (iarg < nargs)
		    xstart = xstart + arglist[iarg+1]
	    }
	}

	return (status)
end

# MAKEBOX -- Compute the coordinates of the corners of a rectangle

procedure makebox (nix, niy, xstart, ystart, xlen, ylen, box)

int	nix		# i: first dimension of output buffer
int	niy		# i: second dimension of output buffer 
real	xstart		# i: x coord lower left corner
real	ystart		# i: y coord lower left corner
real	xlen		# i: x length of box
real	ylen		# i: y length of box
real	box[2,5]	# o: coordinates of box corners
#--
real	xcen, ycen

begin
	# Compute center of output buffer

	xcen = 0.5 * (nix + 1)
	ycen = 0.5 * (niy + 1)

	# lower left corner

	box[1,1] = xcen + xstart
	box[2,1] = ycen + ystart

	# upper left corner

	box[1,2] = xcen + xstart
	box[2,2] = ycen + ystart + ylen

	# upper right corner

	box[1,3] = xcen + xstart + xlen
	box[2,3] = ycen + ystart + ylen

	# lower right corner

	box[1,4] = xcen + xstart + xlen
	box[2,4] = ycen + ystart

	# lower left corner again

	box[1,5] = xcen + xstart
	box[2,5] = ycen + ystart
end

# ROTATBOX -- Rotate a the corners of a rectangle through an angle

procedure rotatbox (nix, niy, box, angle)

int	nix		# i: first dimension of output buffer
int	niy		# i: second dimension of output buffer 
real	box[2,5]	# u: coordinates of box corners
real	angle		# i: rotation angle in degrees
#--
int	vert
real	ca, sa, xcen, ycen, xtemp, ytemp

begin
	ca = cos (DEGTORAD (angle))
	sa = sin (DEGTORAD (angle))

	xcen = 0.5 * (nix + 1)
	ycen = 0.5 * (niy + 1)

	do vert = 1, 5 {
	    xtemp = ca * (box[1,vert] - xcen) + sa * (box[2,vert] - ycen)
	    ytemp = -sa * (box[1,vert] - xcen) + ca * (box[2,vert] - ycen)

	    box[1,vert] = xtemp + xcen
	    box[2,vert] = ytemp + ycen
	}

end

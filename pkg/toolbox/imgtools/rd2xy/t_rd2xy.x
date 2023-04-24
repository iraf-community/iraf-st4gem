include	<math.h>

#  rd2xy -- translate RA and Dec to pixel coordiante 
#
#  Description:
#  ------------
#  
#  Input CL parameters:
#  -----------------
#  "infile"		Input science data file template name
#  "ra"			right ascension
#  "dec"		declination
#  "hour"		the input ra is in hours?
#
#  Output CL parameters:
#  -----------------
#  "x"			X pixel coordinate
#  "y"			Y pixel coordinate
#
#  Date		Author			Description
#  ----		------			-----------
#  12-Aug-1991  J.-C. Hsu		design and coding
#------------------------------------------------------------------------------

procedure t_rd2xy ()

pointer	ipin			# input file pointer
char	infile[SZ_FNAME]
double	x, y
bool	hour
	
double	cd[2, 2], cdinv[2, 2]
double	crpix[2], ra0, dec0, xi, eta, det, bottom

double	ra, dec

pointer	immap()
double	imgetd(), clgetd()
bool	clgetb()
#==============================================================================
begin

	# read input CL parameters
	call clgstr ("infile", infile, SZ_FNAME)
	ra = clgetd ("ra")
	dec = clgetd ("dec")
	hour = clgetb ("hour")

	# open input files
	ipin = immap (infile, READ_ONLY, 0)		

	crpix[1] = imgetd (ipin, "CRPIX1")
	crpix[2] = imgetd (ipin, "CRPIX2")
	ra0 = imgetd (ipin, "CRVAL1")
	dec0 = imgetd (ipin, "CRVAL2")
	cd[1, 1] = imgetd (ipin, "CD1_1")
	cd[1, 2] = imgetd (ipin, "CD1_2")
	cd[2, 1] = imgetd (ipin, "CD2_1")
	cd[2, 2] = imgetd (ipin, "CD2_2")

	det = cd[1,1]*cd[2,2] - cd[1,2]*cd[2,1]
	if (det == 0.d0) 
	    call error (1, "singular CD matrix")

	cdinv[1,1] = cd[2,2] / det
	cdinv[1,2] = -cd[1,2] / det
	cdinv[2,1] = -cd[2,1] / det
	cdinv[2,2] = cd[1,1] / det

	# translate (ra, dec) to (x, y)

	ra0 = DEGTORAD(ra0)
	dec0 = DEGTORAD(dec0)
	if (hour) 
	    ra = ra * 15.d0
	ra = DEGTORAD(ra)
	dec = DEGTORAD(dec)

	bottom = sin(dec)*sin(dec0) + cos(dec)*cos(dec0)*cos(ra-ra0)
	if (bottom == 0.d0) 
	    call error (1, "Unreasonable RA/Dec range")

	xi = cos(dec) * sin(ra-ra0) / bottom
	eta = (sin(dec)*cos(dec0) - cos(dec)*sin(dec0)*cos(ra-ra0)) / bottom
	xi = RADTODEG(xi)
	eta = RADTODEG(eta)

	x = cdinv[1, 1] * xi + cdinv[1, 2] * eta + crpix[1]
	y = cdinv[2, 1] * xi + cdinv[2, 2] * eta + crpix[2]

	call printf ("X = %8.2f, Y = %8.2f\n")
	    call pargd (x)
	    call pargd (y)

	# write to output CL parameters
	call clputd ("x", x)
	call clputd ("y", y)

	# close the input file
	call imunmap (ipin)
end

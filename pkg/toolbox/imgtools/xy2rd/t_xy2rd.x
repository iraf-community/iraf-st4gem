include	<math.h>

#  xy2rd -- translate pixel coordiante to RA and Dec
#
#  "The United States -- bounded on the north by the Aurora Borealis, on the 
#  south by the precession of the equinoxes, on the east by the primeval chaos, 
#  and the west by the Day of Judgement."
#							- John Fiske
#
#  Description:
#  ------------
#  The equations used here is slightly different from but equivalent to 
#  those used in xyeq (gasp packeag).  The equation used here will work 
#  if the reference pixel is at either of the poles.
#  
#  Input CL parameters:
#  -----------------
#  "infile"		Input science data file template name
#  "x"			X pixel coordinate
#  "y"			Y pixel coordinate
#  "hms"		output format switch
#
#  Output CL parameters:
#  -----------------
#  "ra"			right ascension
#  "dec"		declination
#
#  Date		Author			Description
#  ----		------			-----------
#  29-Jul-1991  J.-C. Hsu		design and coding
#------------------------------------------------------------------------------

procedure t_xy2rd ()

pointer	ipin			# input file pointer
char	infile[SZ_FNAME]
double	x, y
bool	hms
	
double	cd[2, 2]
double	crpix[2], ra0, dec0, xi, eta
char	rastr[SZ_LINE], decstr[SZ_LINE]

double	ra, dec

pointer	immap()
double	imgetd(), clgetd()
bool	clgetb()
#==============================================================================
begin

	# read input CL parameters
	call clgstr ("infile", infile, SZ_FNAME)
	x = clgetd ("x")
	y = clgetd ("y")
	hms = clgetb ("hms")

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

	# translate (x,y) to (ra, dec)
	xi = cd[1, 1] * (x - crpix[1]) + cd[1, 2] * (y - crpix[2])
	eta = cd[2, 1] * (x - crpix[1]) + cd[2, 2] * (y - crpix[2])

	xi = DEGTORAD(xi)
	eta = DEGTORAD(eta)
	ra0 = DEGTORAD(ra0)
	dec0 = DEGTORAD(dec0)

	ra = atan2(xi, cos(dec0)-eta*sin(dec0)) + ra0
	dec = atan2(eta*cos(dec0)+sin(dec0), 
			sqrt((cos(dec0)-eta*sin(dec0))**2 + xi**2))

	ra = RADTODEG(ra)
	dec = RADTODEG(dec)
	ra = mod (ra, 360.d0)
	if (ra < 0.d0)
	    ra = ra + 360.d0


	if (hms) {
	    call printf ("RA = %0.4h, Dec = %0.3h\n")
		call pargd (ra/15.d0)
		call pargd (dec)
	    call sprintf (rastr, SZ_LINE, "%0.4h")
		call pargd (ra/15.d0)
	    call sprintf (decstr, SZ_LINE, "%0.3h")
		call pargd (dec)
	} else {
	    call printf ("RA = %11.7f, Dec = %11.7f\n")
		call pargd (ra)
		call pargd (dec)
	    call sprintf (rastr, SZ_LINE, "%11.7f")
		call pargd (ra)
	    call sprintf (decstr, SZ_LINE, "%11.7f")
		call pargd (dec)
	}

	# write to output CL parameters
	call clpstr ("ra", rastr)
	call clpstr ("dec", decstr)

	# close the input file
	call imunmap (ipin)
end

include	<math.h>

#* HISTORY *
#* B.Simon	17-Feb-95	original

# COORDTRANS -- Compute detector coordinate transformation

procedure coordtrans (skycoord, det_ra, det_dec, det_ang, 
		      apscale, apx, apy, mw)

bool	skycoord	# i: use sky coordinates (ra & dec) ?
double	det_ra		# i: right ascension of aperture center (degrees)
double	det_dec		# i: declination of aperture center (degrees)
double	det_ang		# i: position angle of aperture (degrees)
double	apscale		# i: aperture scale (degrees per pixel)
int	apx		# i: aperture x length
int	apy		# i: aperture y length
pointer	mw		# o: world coordinate descriptor
#--
int	ndim
double	pa, r[2], w[2], cd[2,2]
pointer	bufptr

pointer	mw_open()

begin
	# Compute pixel coordinates of reference point

	r[1] = 0.5 * (apx + 1)
	r[2] = 0.5 * (apy + 1)

	# Compute world coordinates of reference point
	# If sky coordinates are not used, 
	# set detector center to origin

	if (skycoord) {
	    pa = DEGTORAD(det_ang)
	    w[1] = det_ra
	    w[2] = det_dec
	} else {
	    pa = 0.0
	    w[1] = 0.0
	    w[2] = 0.0
	}

	# Compute cd matrix

	cd[1,1] = apscale * cos (pa)
	cd[2,1] = - apscale * sin (pa)
	cd[1,2] = apscale * sin (pa)
	cd[2,2] = apscale * cos (pa)

	# RA axis is inverted from usual convention

	if (skycoord) {
	    cd[1,1] = - cd[1,1]
	    cd[2,1] = - cd[2,1]
	}
	    
	# Set up coordinate transformation

	ndim = 2
	bufptr = NULL
	mw = mw_open (bufptr, ndim)

	call mw_newsystem (mw, "world", ndim)

	call mw_swtermd (mw, r, w, cd, ndim)
	call mw_swattrs (mw, 1, "axtype", "ra")
	call mw_swattrs (mw, 1, "wtype", "tan")
	call mw_swattrs (mw, 2, "axtype", "dec")
	call mw_swattrs (mw, 2, "wtype", "tan")

end

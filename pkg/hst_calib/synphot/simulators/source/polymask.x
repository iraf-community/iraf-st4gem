# POLYMASK -- Convert a polygon description to a mask

procedure polymask (npoly, poly, ox, oy, nsub, ntx, nty, mask, flag)

int	npoly		# i: number of polygon vertices
real	poly[2,npoly]	# i: convex polygon vertices
real	ox		# i: x coord of mask center
real	oy		# i: y coord of mask center
int	nsub		# i: number of subpixels per pixel
int	ntx		# i: number x pixels in mask
int	nty		# i: number y pixels in mask
real	mask[ntx,nty]	# u: polygon mask
int	flag		# o: flag indicating if polygon is on mask
#--
int	ity, itx, isub, ip, nsect, ixl, ixr
pointer	sp, sect
real	mx0, my0, mdx, mdy, px0, py0, pdx, pdy
real	step, det, dx0, dy0, dinv, t, s, left, right, frac

string	intercept  "polymask: bad number of intercepts btw mask and polygon"

begin
	# Allocate memory for temporary array

	call smark (sp)
	call salloc (sect, npoly, TY_REAL)

	# Compute horizontal line between edges of mask

	mx0 = int (ox + 0.5) - ntx / 2
	my0 = int (oy + 0.5) - nty / 2
	mdx = ntx
	mdy = 0.0
	step = 1.0 / real(nsub)

	flag = NO

	do ity = 1, nty {
	    do isub = 1, nsub {
		# Find the intersection of the polygon with this line

		nsect = 0
		do ip = 1, npoly-1 {
		    # Calulate line for polygon side

		    px0 = poly[1,ip]
		    py0 = poly[2,ip]
		    pdx = poly[1,ip+1] - poly[1,ip]
		    pdy = poly[2,ip+1] - poly[2,ip]

		    det = pdx * mdy - mdx * pdy
		    if (abs (det) > 1.0e-20) {
			dx0 = px0 - mx0
			dy0 = py0 - my0
			dinv = 1.0 / det

			# If t is btw 0 and 1, the line intersects 
			# the polygon side

			t = dinv * (mdx * dy0 - mdy * dx0)
			if (t >= 0.0 && t <= 1.0) {
			    s = dinv * (pdx * dy0 - pdy * dx0)
			    Memr[sect+nsect] = s
			    nsect = nsect + 1
			}
		    }
		}

		# If the horizontal line intersects the polygon, set
		# the mask pixels between the intersection

		if (nsect == 2) {
		    left = min (Memr[sect], Memr[sect+1])
		    right = max (Memr[sect], Memr[sect+1])

		    if (right >= 0.0 && left <= 1.0) {
			flag = YES
			left = max (left, 0.0)
			right = min (right, 1.0)

			ixl = left * ntx + 1.9999
			ixr = right * ntx

			# First, set pixels entirely within the mask

			do itx = ixl, ixr {
			    mask[itx,ity] = mask[itx,ity] + step
			    mask[itx,ity] = min (1.0, mask[itx,ity])
			}

			# Handle endpoint pixels separately

			if (ixl > 1) {
			    itx = ixl - 1
			    frac =  aint ((itx - left * ntx) * nsub) * step
			    mask[itx,ity] = mask[itx,ity] + step * frac
			    mask[itx,ity] = min (1.0, mask[itx,ity])
			}

			if (ixr < ntx) {
			    itx = ixr + 1
			    frac =  aint ((right * ntx - ixr) * nsub) * step
			    mask[itx,ity] = mask[itx,ity] + step * frac
			    mask[itx,ity] = min (1.0, mask[itx,ity])
			}
		    }

		} else if (nsect != 0) {
		    # Check for algorithm error
		    call printerr_int (intercept, nsect)
		}

		# Step horizontal line up to next line of pixels

		my0 = my0 + step
	    }
	}

	call sfree (sp)
end

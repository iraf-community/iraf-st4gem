#* HISTORY *
#* B.Simon	06-Nov-95	original

# PSFSHIFT -- Shift the psf to account for object placement within pixel

procedure psfshift (ox, oy, nsub, npx, npy, in, out)

real	ox		# i: object x position
real	oy		# i: object y position
int	nsub		# i: number of subpixels per pixel
int	npx		# i: x dimension of psf arrays
int	npy		# i: y dimension of psf arrays
real	in[npx,npy]	# i: input array
real	out[npx,npy]	# o: output array
#--
int	ix, iy, jx, jy, kx, ky
real	xpos, ypos, xshift, yshift, xfac[3], yfac[3]

begin
	# Compute interpolation factors from offset in object position

	xpos = nsub * (ox + 0.5 - int (ox + 0.5))
	xshift = (xpos - int (xpos)) - 0.5
	xfac[1] = abs (min (xshift, 0.0))
	xfac[3] = max (xshift, 0.0)
	xfac[2] = 1.0 - xfac[1] - xfac[3]

	ypos = nsub * (oy + 0.5 - int (oy + 0.5))
	yshift = (ypos - int (ypos)) - 0.5
	yfac[1] = abs (min (yshift, 0.0))
	yfac[3] = max (yshift, 0.0)
	yfac[2] = 1.0 - yfac[1] - yfac[3]

	# Multiply psf by interpolation factors to get shifted psf

	do iy = 1, npy {
	    do ix = 1, npx {
		out[ix,iy] = 0.0
		do jy = -1, 1 {
		    ky = iy + jy
		    if (ky < 1 || ky > npy || yfac[jy+2] <= 0.0)
			next

		    do jx = -1, 1 {
			kx = ix + jx
			if (kx < 1 || kx > npx || xfac[jx+2] <= 0.0)
			    next

			out[ix,iy] = out[ix,iy] + xfac[jx+2] * yfac[jy+2] *
				     in[kx,ky]
		    }
		}
	    }
	}

end


#* HISTORY *
#* B.Simon	17-Feb-95	original
#* B.Simon	04-Feb-98	replaced oob error with next

# STARCONV -- Produce an image of a star convolved with a psf

procedure starconv (ox, oy, flux, nsub, npx, npy, ntx, nty, apsf, anobj)

real	ox		# i: star x position
real	oy		# i: star y position
real	flux		# i: star flux
int	nsub		# i: number of pixel subdivisions
int	npx		# i: x dimension of psf
int	npy		# i: y dimension of psf
int	ntx		# i: x dimension of convolved star
int	nty		# i: y dimension of convolved star
real	apsf[npx,npy]	# i: point spread function
real	anobj[ntx,nty]	# o: convolved star
#--
int	jsx, jsy, itx, ity, isx, isy, ipx, ipy

string	convolve "Convolving with PSF"

begin
	# Set convolved star array to zero

	call aclrr (anobj, ntx*nty)

	# Compute offset to star location within pixel

	jsx = (ox + 0.5 - int (ox + 0.5)) * nsub + 1
	jsy = (oy + 0.5 - int (oy + 0.5)) * nsub + 1

	# Compute location of far edge of the convolved star

	jsx = nsub * (ntx / 2) + (npx + 1) / 2 + jsx - 1
	jsy = nsub * (nty / 2) + (npy + 1) / 2 + jsy - 1

	# Break location into the pixel number and the subpixel number

	ity = jsy / nsub
	isy = jsy - nsub * ity
	if (isy == 0) {
	    isy = nsub
	} else {
	    ity = ity + 1
	}

	do ipy = 1, npy {
	    if (ity > nty || ity < 1)
		next

	    # Print diagnostic message

	    call done_message (convolve, ipy-1, npy)

	    # Break location into the pixel number and the subpixel number

	    itx = jsx / nsub
	    isx = jsx - nsub * itx
	    if (isx == 0) {
		isx = nsub
	    } else {
		itx = itx + 1
	    }

	    # Sum the psf flux lying within the current convolved pixel

	    do ipx = 1, npx {
		if (itx > ntx || itx < 1)
		    next

		anobj[itx,ity] = anobj[itx,ity] + apsf[ipx,ipy]

		if (isx > 1) {
		    isx = isx - 1
		} else {
		    isx = nsub
		    itx = itx - 1
		}
	    }

	    if (isy > 1) {
		isy = isy - 1
	    } else {
		isy = nsub
		ity = ity - 1
	    }
	}

	# Multiply the (normalized) star image by the flux

	call amulkr (anobj, flux, anobj, ntx*nty)

	# Finish diagnostic message

	call done_message (convolve, npy, npy)

end

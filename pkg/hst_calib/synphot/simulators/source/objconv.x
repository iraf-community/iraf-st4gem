#* HISTORY *
#* B.Simon	17-Feb-95	original
#* B.Simon	04-May-95	Changed shape descriptor
#* B.Simon	04-Feb-98	replaced oob error with next

# OBJCONV -- Convolve an extended object with a psf

procedure objconv (ox, oy, flux, nsub, nox, noy, npx, npy, ntx, nty, 
		   valfunc, apsf, anobj)

real	ox		# i: object x center
real	oy		# i: object y center
real	flux		# i: object flux
int	nsub		# i: number of pixel subdivisions
int	nox		# i: x dimension of shape
int	noy		# i: y diemsion of shape
int	npx		# i: x dimension of psf
int	npy		# i: y dimension of psf
int	ntx		# i: x dimension of convolved object
int	nty		# i: y dimension of convolved object
pointer	valfunc		# i: function to compute shape value
real	apsf[npx,npy]	# i: point spread function
real	anobj[ntx,nty]	# o: convolved object
#--
int	jsx, jsy, itx, ity, isx, isy, iox, ioy, ipx, ipy
real	ashp, sum, total

string	convolve "Convolving with PSF"

begin
	# Set convolved star array to zero

	call aclrr (anobj, ntx*nty)

	# Compute offset to star location within pixel

	jsx = (ox + 0.5 - int (ox + 0.5)) * nsub + 1
	jsy = (oy + 0.5 - int (oy + 0.5)) * nsub + 1

	# Compute location of near edge of the convolved star

	jsx = nsub * (ntx / 2) - (nox + npx + 1) / 2 + jsx - 1
	jsy = nsub * (nty / 2) - (noy + npy + 1) / 2 + jsy - 1

	sum = INDEFR
	total = 0.0

	do ioy = 1, noy {
	    # Print diagnostic message
	    call done_message (convolve, ioy-1, noy)

	    do iox = 1, nox {
		# Compute shape value at this pixel

		call zcall4 (valfunc, iox, ioy, sum, ashp)
		total = total + ashp

		# Compute location of top edge of convolution

		isy = ioy + npy + jsy

		ity = isy / nsub
		isy = isy - nsub * ity
		if (isy == 0) {
		    isy = nsub
		} else {
		    ity = ity + 1
		}

		# Compute the flux lying within the current convolved pixel

		do ipy = 1, npy  {
		    if (ity > nty || ity < 1)
			next

		    # Compute location of left edge of convolution

		    isx = iox + npx + jsx

		    itx = isx / nsub 
		    isx = isx - nsub * itx
		    if (isx == 0) {
			isx = nsub
		    } else {
			itx = itx + 1
		    }

		    do ipx = 1, npx {
			if (itx > ntx || itx < 1)
			    next

			anobj[itx,ity] = anobj[itx,ity] + apsf[ipx,ipy] *
					 ashp

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
	    }
	}

	# Multiply the (normalized) object image by the flux

	if (total <= 0.0)
	    total = 1.0

	call amulkr (anobj, flux / total, anobj, ntx*nty)

	# Finish diagnostic message

	call done_message (convolve, noy, noy)

end

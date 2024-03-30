#* HISTORY *
#* B.Simon	13-Jul-95	original
#* B.Simon	04-Feb-98	replaced oob error with next

# LSFCONV -- Convolve a source with a line spread function

procedure lsfconv (flux, nsub, nlsf, alsf, ntx, nty, anobj, ndx, ndy, asrc)

real	flux		# i: integrated flux of source
int	nsub		# i: number of subpixels per pixel
int	nlsf		# i: length of lsf
real	alsf[nlsf]	# i: line spread function (lsf)
int	ntx		# i: x dimension of unconvolved source
int	nty		# i: y dimension of unconvolved source
real	anobj[ntx,nty]	# i: unconvolved source
int	ndx		# i: x dimension of convolved source
int	ndy		# i: y dimension of convolved source
real	asrc[ndx,ndy]	# o: convolved source
#--
int	itx, ity, isy, jsy, joy, ilsf
real	total

real	asumr()

begin
	# Set output array to zero

	call aclrr (asrc, ndx*ndy)

	do itx = 1, ntx {
	    do ity = 1, nty {
		# Compute offset to initial pixel in object

		isy = nlsf + nsub * (ity - 1)
		joy = isy / nsub
		jsy = isy - nsub * joy
		if (jsy == 0) {
		    jsy = nsub
		} else {
		    joy = joy + 1
		}

		# Convolve lsf with a line of the object

		do ilsf = 1, nlsf {
		    if (joy > ndy || joy < 1)
			next

		    asrc[itx,joy] = asrc[itx,joy] + alsf[ilsf] * anobj[itx,ity]

		    jsy = jsy - 1
		    if (jsy == 0) {
			joy = joy - 1
			jsy = nsub
		    }
		}
	    }
	}

	# Normalize convolved image to flux

	total = asumr (asrc, ndx*ndy)
	if (total <= 0.0)
	    total = 1.0

	call amulkr (asrc, flux / total, asrc, ndx*ndy)

end

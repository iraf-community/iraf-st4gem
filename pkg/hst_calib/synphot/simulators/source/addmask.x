# ADDMASK -- Combine source with aperture mask

procedure addmask (ntx, nty, mask, asrc, flux)

int	ntx		# i: first dimension of mask and source
int	nty		# i: second dimension of mask and source
real	mask[ntx,nty]	# i: aperture mask
real	asrc[ntx,nty]	# u: source
real	flux		# o: flux of masked source
#--
int	itx, ity

begin
	# Multiply source by apeture mask

	flux = 0.0
	do ity = 1, nty {
	    do itx = 1, ntx {
		asrc[itx,ity] = asrc[itx,ity] * mask[itx,ity]
		flux = flux + asrc[itx,ity]
	    }
	}

end

#* HISTORY *
#* B.Simon	17-Feb-95	original

# ADDOBJECT -- Add an object to the output array

procedure addobject (ox, oy, ntx, nty, nix, niy, anobj, out)

real	ox		# i: object x center
real	oy		# i: object y center
int	ntx		# i: x dimension of object
int	nty		# i: y dimension of object
int	nix		# i: x dimension of output
int	niy		# i: y dimension of output
real	anobj[ntx,nty]	# i: convolved object
real	out[nix,niy]	# o: output array
#--
int	iix, iiy, jix, jiy, kix, kiy, itx, ity, jtx, jty

begin
	# Clip object to output array boundaries

	jtx = 1
	jix = int (ox + 0.5) - ntx / 2
	kix = min (jix + (ntx - 1), nix)

	if (jix < 1) {
	    jtx = 1 + (1 - jix)
	    jix = 1
	}
	
	jty = 1
	jiy = int (oy + 0.5) - nty / 2
	kiy = min (jiy + (nty - 1), niy)

	if (jiy < 1) {
	    jty = 1 + (1 - jiy)
	    jiy = 1
	}

	# Add object to output array 

	ity = jty
	do iiy = jiy, kiy  {
	    itx = jtx

	    do iix = jix, kix {
		out[iix,iiy] = out[iix,iiy] + anobj[itx,ity]
		itx = itx + 1
	    }

	    ity = ity + 1
	}

end

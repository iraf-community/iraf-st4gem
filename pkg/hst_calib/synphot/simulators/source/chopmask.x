#* HISTORY *
#* B.Simon	18-Jul-95	original

# CHOPMASK -- Replace mask with smaller one containing region of interest

procedure chopmask (ntx, nty, mask)

int	ntx		# u: x dimension of mask
int	nty		# u: y dimension of mask
pointer	mask		# u: pointer to mask array
#--
int	itx, ity, ixlo, ixhi, iylo, iyhi, xdim, ydim
pointer	ptr, ptr1, ptr2, mask2

string	nomask  "Aperture mask not found"

begin
	# Find box containing nonzero portion of mask

	ixlo = ntx + 1
	ixhi = 0
	iylo = nty + 1
	iyhi = 0
	ptr = mask

	do ity = 1, nty {
	    do itx = 1, ntx {
		if (Memr[ptr] > 0.0) {
		    ixlo = min (ixlo, itx)
		    ixhi = max (ixhi, itx)
		    iylo = min (iylo, ity)
		    iyhi = max (iyhi, ity)
		}

		ptr = ptr + 1
	    }
	}

	if (ixhi == 0)
	    call printerr_str (nomask, "chopmask")

	# Allocate new mask and copy boxed region into it

	xdim = ixhi - ixlo + 1
	ydim = iyhi - iylo + 1

	call malloc (mask2, xdim*ydim, TY_REAL)

	ptr2 = mask2
	ptr1 = mask + ntx * (iylo - 1)

	do ity = iylo, iyhi {
	    ptr = ptr1 + ixlo - 1

	    do itx = ixlo, ixhi {
		Memr[ptr2] = Memr[ptr]
		ptr2 = ptr2 + 1
		ptr = ptr + 1
	    }

	    ptr1 = ptr1 + ntx
	}

	# Free old mask and overwrite subroutine parameters

	call mfree (mask, TY_REAL)
	mask = mask2
	ntx = xdim
	nty = ydim

end

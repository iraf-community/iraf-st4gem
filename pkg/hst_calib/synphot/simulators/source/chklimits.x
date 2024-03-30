include	<mach.h>
include	"limit.h"

#* HISTORY *
#* B.Simon	15-Aug-95	original

# CHKLIMITS -- Check count rate limits and print a warning message if exceeded

procedure chklimits (limit, exptime, nix, niy, out)

real	limit[LEN_L]	# i: count rate limits
real	exptime		# i: exposure time
int	nix		# i: x dimension of output image
int	niy		# i: y dimension of output image
real	out[nix,niy]	# i: output array
#--
int	ix, iy, ncount, nbright
pointer	sp, number, value
real	total, maxcount, maxbright

string	badlcount  "Local count rate limit exceeded. # times"
string	badlbright "Local bright rate limit exceeded. # times"
string	badgcount  "Global count rate limit exceeded. Total flux"
string	badgbright "Global bright rate limit exceeded. Total flux"

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (number, SZ_FNAME, TY_CHAR)
	call salloc (value, SZ_FNAME, TY_CHAR)

	# Check local limits

	if (IS_INDEFR(limit[LCOUNT])) {
	    maxcount = MAX_REAL
	} else {
	    maxcount = exptime * limit[LCOUNT]
	}

	if (IS_INDEFR(limit[LBRIGHT])) {
	    maxbright = MAX_REAL
	} else {
	    maxbright = exptime * limit[LBRIGHT]
	}

	ncount = 0
	nbright = 0
	total = 0.0

	do iy = 1, niy {
	    do ix = 1, nix {
		total = total + out[ix,iy]

		if (out[ix,iy] > maxcount)
		    ncount = ncount + 1

		if (out[ix,iy] > maxbright)
		    nbright = nbright + 1
	    }
	}

	# Print message if local limits exceeded

	if (ncount > 0) {
	    call sprintf (Memc[number], SZ_FNAME, "%d")
	    call pargi (ncount)

	    call synphotwarn (badlcount, Memc[number])
	}

	if (nbright > 0) {
	    call sprintf (Memc[number], SZ_FNAME, "%d")
	    call pargi (nbright)

	    call synphotwarn (badlbright, Memc[number])
	}

	# Check global and bright limits

	total = total / exptime

	if (! IS_INDEFR(limit[GCOUNT])) {
	    if (total > limit[GCOUNT]) {
		call sprintf (Memc[value], SZ_FNAME, "%g")
		call pargr (total)
		call synphotwarn (badgcount, Memc[value])
	    }
	}

	if (! IS_INDEFR(limit[GBRIGHT])) {
	    if (total > limit[GBRIGHT]) {
		call sprintf (Memc[value], SZ_FNAME, "%g")
		call pargr (total)
		call synphotwarn (badgbright, Memc[value])
	    }
	}

	call sfree (sp)
end

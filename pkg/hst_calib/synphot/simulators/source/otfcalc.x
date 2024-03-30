include "otf.h"

#* HISTORY *
#* B.Simon	17-Feb-95	original
#* B.Simon	29-Jun-95	revised so nxpix != nypix

# OTFCALC -- Calculate the weighted OTF function

procedure otfcalc (otf, ox, oy, nsub, nweight, weight, nxpix, nypix, array)

pointer	otf		    # i: otf descriptor
real	ox		    # i: object x position
real	oy		    # i: object y position
int	nsub		    # i: number of subpixels per pixel
int	nweight		    # i: number of weights
real	weight[ARB]	    # i: weighting factors
int	nxpix		    # i: x dimension of otf array
int	nypix		    # i: y dimension of otf array
real	array[nxpix,nypix]  # o: combined otf array
#--
int	idx

string	badweight  "otfcalc: bad value for nweight"
string	badnpix    "otfcalc: bad value for npix"

begin
	# Check inputs and structure for consistency

	if (OTF_NUMBER(otf) != nweight) 
	    call printerr_int (badweight, nweight)

	if (OTF_NXPIX(otf) != nxpix)
	    call printerr_int (badnpix, nxpix)

	if (OTF_NYPIX(otf) != nypix)
	    call printerr_int (badnpix, nypix)

	# Initialize sum

	call amovkr (0.0, array, nxpix*nypix)

	# Sum otfs multiplied by weights

	do idx = 1, nweight {
	    if (weight[idx] <= 0.0)
		next

	    call wtsum (nxpix, nypix, weight[idx], 
			OTF_ARRAY(otf,idx), array)
	}


end

# WTSUM -- Compute weighted sum

procedure wtsum (nxpix, nypix, weight, in, out)

int	nxpix		# i: x dimension of otf arrays
int	nypix		# i: y dimension of otf arrays
real	weight		# i: weighting factor
real	in[nxpix,nxpix]	# i: input array
real	out[nxpix,nypix]# u: output array
#--
int	ix, iy

begin
	do iy = 1, nypix
	    do ix = 1, nxpix
		out[ix,iy] = out[ix,iy] + in[ix,iy] * weight

end

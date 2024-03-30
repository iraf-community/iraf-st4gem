#* HISTORY *
#* B.Simon	08-Jun-94	adopted from hgline
	
# HISTPLOT -- Plot a stepped curve of flux vs. wavelength

procedure histplot (gp, wave, flux, nwave)

pointer	gp		# Graphics descriptor
real	wave[ARB]	# X coordinates data
real	flux[ARB]	# Y coordinates of data
int	nwave		# Number of datapoints
#--
int	iwave
real	x, y

begin
	# Do the first horizontal line

	x = wave[1]
	y = flux[1]
	call gamove (gp, x, y)
	x = wave[2]
	call gadraw( gp, x, y )

	do iwave = 2, nwave-1 {
	    y = flux[iwave]

	    # vertical connection
	    call gadraw (gp, x, y)

	    # horizontal line
	    x = wave[iwave+1]
	    call gadraw (gp, x, y)
	}

	# Draw last segments

	y = flux[nwave]
	call gadraw (gp, x, y)
end

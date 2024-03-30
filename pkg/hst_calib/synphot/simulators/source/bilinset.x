#* HISTORY *
#* B.Simon	17-Feb-95	original

# BILINSET -- Set the parameters used by the blinear interpolator

procedure bilinset (buffer, xlen, ylen)

pointer	buffer		# i: Buffer to be inerpolated on
int	xlen		# i: First dimension of buffer
int	ylen		# i: Second dimension of buffer
#--
include "bilin.com"

begin
	buf = buffer
	nx = xlen
	ny = ylen
end


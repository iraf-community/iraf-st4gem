#* HISTORY *
#* B.Simon	11-Jul-95	derived from bilinset

# LINEARSET -- Set the parameters used by the linear interpolator

procedure linearset (buffer, xlen)

pointer	buffer		# i: Buffer to be inerpolated on
int	xlen		# i: Length of buffer
#--
include "linear.com"

begin
	buf = buffer
	nx = xlen
end


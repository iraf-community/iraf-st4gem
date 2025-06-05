include <iraf77.h>

# UHDCKL -- Close the image header field name list and return all
# associated storage

procedure uhdckl (fn, istat)

pointer	fn
int	istat		# return error code

begin
	istat = ER_OK
	call imcfnl (fn)
end

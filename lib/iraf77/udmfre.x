include <iraf77.h>

# UDMFRE -- free dynamic memory from the VOS

procedure udmfre (ptr, dtype, status)

pointer ptr		#buffer pointer 
int	dtype		#datatype of the storage elements
int	status		#return status

begin
	status = ER_OK
	iferr (call mfree (ptr, dtype))
		status = ER_FREEDYNMEM

end

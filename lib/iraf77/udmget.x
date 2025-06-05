include <iraf77.h>

# UDMGET -- get dynamic memory from the VOS

procedure udmget (nelem, dtype, ptr, status)

int	nelem		#number of elements of storage required
int	dtype		#datatype of the storage elements
pointer ptr		#buffer pointer for output
int	status		#return status

# Question - what datatypes do we support here? 
# IRAF char type will cause problems in FORTRAN programs.

begin
	status = ER_OK

	if (dtype == TY_CHAR)
		status = ER_DMINVDTYP

	else   #if (dtype == other types?)
		iferr (call malloc (ptr, nelem, dtype))
			status = ER_DYNMEMALC

end

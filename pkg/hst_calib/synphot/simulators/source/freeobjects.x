include "object.h"

#* HISTORY *
#* B.Simon	17-Feb-95	original

# FREE_OBJECTS -- Free memory allocated for object descriptor

procedure free_objects (obj)

pointer	obj		# i: object descriptor
#--
int	index
pointer	shape

begin
	do index = 1, OBJ_NUMBER(obj) {
	    shape = OBJ_SHAPE(obj,index)
	    if (shape != NULL)
		call free_func (shape)
	}

	call mfree (OBJ_SHPARY(obj), TY_INT)
	call mfree (OBJ_PWTARY(obj), TY_REAL)
	call mfree (OBJ_FLUXARY(obj), TY_REAL)
	call mfree (OBJ_YPOSARY(obj), TY_REAL)
	call mfree (OBJ_XPOSARY(obj), TY_REAL)

	call mfree (obj, TY_INT)
end

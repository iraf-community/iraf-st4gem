include "function.h"

#* HISTORY *
#* B.Simon	17-Feb-95	original
#* B.Simon	04-May-95	Changed function descriptor

# FREE_FUNC -- Free memory allocated to  function strcucture

procedure free_func (func)

pointer	func		# i: function descriptor
#--
int	ipar

begin
	# Check to see if function was defined

	if (func == NULL)
	    return

	# Free parameters

	do ipar = 1, FUN_NPAR(func)
	    call mfree (FUN_PARAM(func,ipar), FUN_TYPE(func,ipar))

	# Free parameter and type arrays

	call mfree (FUN_PARARY(func), TY_INT)
	call mfree (FUN_TYPARY(func), TY_INT)

	# Free function structure

	call mfree (func, TY_INT)
end

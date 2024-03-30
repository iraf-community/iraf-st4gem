include "dtrans.h"

#* HISTORY *
#* B.Simon	22-Jul-94	original

# CLSDTRANS -- Release the spectrum descriptor

procedure clsdtrans (dtrans)

pointer	dtrans		# u: spectrum descriptor
#--

begin
	if (dtrans == NULL)
	    return

	call mfree (TRN_XEFF(dtrans), TY_REAL)
	call mfree (TRN_YEFF(dtrans), TY_REAL)

	call mfree (dtrans, TY_STRUCT)
	dtrans = NULL
end

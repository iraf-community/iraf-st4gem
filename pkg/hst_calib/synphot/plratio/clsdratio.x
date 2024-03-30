include "dratio.h"

#* HISTORY *
#* B.Simon	12-Jul-94	original

# CLSDRATIO -- Release the spectrum descriptor

procedure clsdratio (dratio)

pointer	dratio		# u: spectrum descriptor
#--

begin
	if (dratio == NULL)
	    return

	call mfree (RAT_WAVE(dratio), TY_REAL)
	call mfree (RAT_FLUX(dratio), TY_REAL)

	call mfree (dratio, TY_STRUCT)
	dratio = NULL
end

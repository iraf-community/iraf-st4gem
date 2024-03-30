include "dspec.h"

#* HISTORY *
#* B.Simon	13-Jun-94	original

# CLSDSPEC -- Release the spectrum descriptor

procedure clsdspec (dspec)

pointer	dspec		# u: spectrum descriptor
#--

begin
	if (dspec == NULL)
	    return

	call mfree (SPC_WAVE(dspec), TY_REAL)
	call mfree (SPC_FLUX(dspec), TY_REAL)
	call mfree (SPC_STIM(dspec), TY_REAL)
	call mfree (SPC_PIVOT(dspec), TY_REAL)
	call mfree (SPC_FWHM(dspec), TY_REAL)

	call mfree (dspec, TY_STRUCT)
	dspec = NULL
end

include "dband.h"

#* HISTORY *
#* B.Simon	10-Jun-94	original

# CLSDBAND -- Release the bandpass descriptor

procedure clsdband (dband)

pointer	dband		# u: bandpass descriptor
#--

begin
	if (dband == NULL)
	    return

	call mfree (BND_BAND(dband), TY_REAL)
	call mfree (BND_WAVE(dband), TY_REAL)

	call mfree (dband, TY_STRUCT)
	dband = NULL
end

include "dphot.h"

#* HISTORY *
#* B.Simon	15-Jun-94	original

# CLSDPHOT -- Release the photometry descriptor

procedure clsdphot (dphot)

pointer	dphot		# u: spectrum descriptor
#--

begin
	if (dphot == NULL)
	    return

	call mfree (PHT_FWHM(dphot), TY_REAL)
	call mfree (PHT_PIVOT(dphot), TY_REAL)
	call mfree (PHT_STIM(dphot), TY_REAL)
	call mfree (PHT_BAND(dphot), TY_REAL)
	call mfree (PHT_WAVE(dphot), TY_REAL)

	call mfree (dphot, TY_STRUCT)
	dphot = NULL
end

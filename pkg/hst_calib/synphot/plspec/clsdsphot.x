include "dsphot.h"

#* HISTORY *
#* B.Simon	13-Jun-94	original

# CLSDSPHOT -- Release the spectrophotometry descriptor

procedure clsdsphot (dsphot)

pointer	dsphot		# u: spectrophotometry descriptor
#--

begin
	if (dsphot == NULL)
	    return

	call mfree (SPT_FWHM(dsphot), TY_REAL)
	call mfree (SPT_ERROR(dsphot), TY_REAL)
	call mfree (SPT_FLUX(dsphot), TY_REAL)
	call mfree (SPT_WAVE(dsphot), TY_REAL)

	call mfree (dsphot, TY_STRUCT)
	dsphot = NULL
end

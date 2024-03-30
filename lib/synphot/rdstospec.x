include	<tbset.h>
include "libsynphot.h"

# RDSTOSPEC -- Read a stored spectrum, skip error checking

procedure rdstospec (fname, nwave, wave, spec)

char	fname[ARB]	# i: table file name
int     nwave           # i: length of wavelength and spectrum arrays
real    wave[ARB]       # i: wavelengths at which spectrum is computed
real	spec[ARB]	# o: spectrum flux
#--
int	nrow, irow, utype
pointer	sp, wavval, fluxval, nulflg, units
pointer	tp, wv, fx
real	factor

string	wavecol   "WAVELENGTH"
string	fluxcol   "FLUX"
string	defunits  "flam"
string  shortlist "photlam,flam"
string  badunits  "Flux units not supported for standard files"

int	tbpsta(), word_match()
pointer	opnsyntab()

begin
	# Open table containg spectrum

	tp = opnsyntab (fname)

	# Allocate temporary arrays to hold contents of table

	nrow = tbpsta (tp, TBL_NROWS)

	call smark (sp)
	call salloc (wavval, nrow, TY_REAL)
	call salloc (fluxval, nrow, TY_REAL)
	call salloc (nulflg, nrow, TY_BOOL)
	call salloc (units, SZ_COLUNITS, TY_CHAR)

	# Read wavelength array

	call syncolptr (tp, wavecol, 1, wv)
	call tbcgtr (tp, wv, Memr[wavval], Memb[nulflg], 1, nrow)

	# Read flux array

	call syncolptr (tp, fluxcol, 2, fx)
	call tbcgtr (tp, fx, Memr[fluxval], Memb[nulflg], 1, nrow)
	call tbcigt (fx, TBL_COL_UNITS, Memc[units], SZ_COLUNITS)

	# Set default flux units if not found in table

	call strfix (Memc[units])
	if (Memc[units] == EOS)
	    call strcpy (defunits, Memc[units], SZ_COLUNITS)

	# Convert flux units to PHOTLAM. This is done in-line 
	# because calling anytophot causes function recursion

	call strfix (Memc[units])
	utype = word_match (Memc[units], shortlist)
	switch (utype) {
	case 1:  # photlam (no change)
	    ;

	case 2:  # flam
	    factor = 1.0 / (H * C) 
	    do irow = 0, nrow - 1 {
		Memr[fluxval+irow] = factor * Memr[wavval+irow] * 
				     Memr[fluxval+irow]
	    }

	default:
	    call synphoterr (badunits, Memc[units])
	}
	

	# Interpolate on wavelength grid

	call syninterp (nrow, Memr[wavval], Memr[fluxval], 
			nwave, wave, spec)

	call sfree (sp)
end

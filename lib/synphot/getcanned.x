include	<tbset.h>

# GETCANNED -- Read a canned spectrum, skip error checking for speed

procedure getcanned (fname, nwave, wave, spec)

char	fname[ARB]	# i: table file name
int     nwave           # i: length of wavelength and spectrum arrays
real    wave[ARB]       # i: wavelengths at which spectrum is computed
real	spec[ARB]	# o: spectrum flux
#--
int	nrow
pointer	sp, wavval, fluxval, nulflg
pointer	tp, wv, fx
real	extrap

string	wavecol   "WAVELENGTH"
string	fluxcol   "FLUX"

int	tbpsta()
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

	# Read wavelength array

	call syncolptr (tp, wavecol, 1, wv)
	call tbcgtr (tp, wv, Memr[wavval], Memb[nulflg], 1, nrow)

	# Read flux array

	call syncolptr (tp, fluxcol, 2, fx)
	call tbcgtr (tp, fx, Memr[fluxval], Memb[nulflg], 1, nrow)

	# Interpolate on wavelength grid

	extrap = min (Memr[fluxval], Memr[fluxval+nrow-1])
	call synextrap (extrap, nrow, Memr[wavval], Memr[fluxval], 
			nwave, wave, spec)

	call sfree (sp)
end

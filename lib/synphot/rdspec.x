include	<tbset.h>

#* HISTORY *
#* B.Simon	22-Jan-97	syninterp replaced by rebin

# RDSPEC -- Read a spectrum from a table sampled on the wavelength set

procedure rdspec (filename, nwave, wave, spec)

char	filename[ARB]	# i: name of file containing spectrum
int     nwave           # i: length of wavelength and spectrum arrays
real    wave[ARB]       # i: wavelengths at which spectrum is computed
real	spec[ARB]	# o: spectral flux at sampled wavelengths
#--
int	ic, nc, irow, nrow, done, fluxnum
pointer	sp, fluxname, tabname, units, wavval, fluxval, nulflg
pointer	tp, wv, fx
real	fillval

string	wavecol   "WAVELENGTH"
string	fluxcol   "FLUX"
string  defunits  "flam"
string	nullwave  "Indef found in wavelength column"
string	badwave   "Wavelength units not found, angstroms assumed"
string	badflux   "Flux units not found, FLAM assumed"
string	nullflux  "Flux column all indef"
string	wavorder  "Wavelength column is not sorted"

int	ctoi(), tbpsta(), anytoang(), anytophot()
int	fillnull(), is_magunit(), synsort2()
pointer	opnsyntab()

errchk	opnsyntab, syncolptr, anytophot, synphoterr

begin
	# Allocate memory for temporary string

	call smark (sp)
	call salloc (fluxname, SZ_FNAME, TY_CHAR)
	call salloc (tabname, SZ_FNAME, TY_CHAR)
	call salloc (units, SZ_COLUNITS, TY_CHAR)

	# Extract table and column name for file name and open table

	call breakcomp (filename, Memc[tabname], Memc[fluxname], SZ_FNAME)
	tp = opnsyntab (Memc[tabname])

	# Allocate temporary arrays to hold contents of table

	nrow = tbpsta (tp, TBL_NROWS)
	call salloc (wavval, nrow, TY_REAL)
	call salloc (fluxval, nrow, TY_REAL)
	call salloc (nulflg, nrow, TY_BOOL)

	# Read wavelength array

	call syncolptr (tp, wavecol, 1, wv)
	call tbcgtr (tp, wv, Memr[wavval], Memb[nulflg], 1, nrow)

	# Check for nulls

	do irow = 1, nrow {
	    if (Memb[nulflg+irow-1])
		call synphoterr (nullwave, filename)
	}

	# Convert wavelength units to angstroms

	call tbcigt (wv, TBL_COL_UNITS, Memc[units], SZ_COLUNITS)

	if (anytoang (Memc[units], Memr[wavval], nrow) == NO) {
	    if (tbpsta (tp, TBL_WHTYPE) != TBL_TYPE_TEXT)
		call synphotwarn (badwave, filename)
	}

	# Get default flux column name if not found in file name

	if (Memc[fluxname] == EOS) {
	    fluxnum = 2
	    call strcpy (fluxcol, Memc[fluxname], SZ_FNAME)

	} else {
	    ic = 1
	    nc = ctoi (Memc[fluxname], ic, fluxnum)
	}

	# Read flux array and remove nulls

	call syncolptr (tp, Memc[fluxname], fluxnum, fx)
	call tbcgtr (tp, fx, Memr[fluxval], Memb[nulflg], 1, nrow)
	call tbcigt (fx, TBL_COL_UNITS, Memc[units], SZ_COLUNITS)

	if (is_magunit (Memc[units]) == NO) {
	    fillval = 0.0
	} else {
	    fillval = 100.0
	}

	if (fillnull (fillval, nrow, Memr[wavval], Memr[fluxval]) == ERR)
	    call synphoterr (nullflux, filename)

	# Set default flux units if not found in table

	call strfix (Memc[units])
	if (Memc[units] == EOS) {
	    call strcpy (defunits, Memc[units], SZ_COLUNITS)

	    if (tbpsta (tp, TBL_WHTYPE) != TBL_TYPE_TEXT)
		call synphotwarn (badflux, filename)
	}

	# Convert flux units to PHOTLAM

	done = anytophot (Memc[units], nrow, Memr[wavval], Memr[fluxval])

	# Sort flux and wavelength in ascending order 
	# and rebin on wavelength grid

	if (synsort2 (nrow, Memr[wavval], Memr[fluxval]) == ERR)
	    call synphoterr (wavorder, filename)

	call rebin (nrow, Memr[wavval], Memr[fluxval], 
		    nwave, wave, spec)

	call sfree (sp)
end

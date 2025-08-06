include	<tbset.h>

#* HISTORY *
#* B.Simon	10-Jun-94	Original

# RDSPHOT -- Read a spectrophotmetry table

procedure rdsphot (table, nwave, wave, flux, err, fwhm)

char	table[ARB]	# i: spectrophotometry table name
int	nwave		# i: lnegth of wavelength and other arrays
real	wave[ARB]	# i: wavelength array
real	flux[ARB]	# o: flux array
real	err[ARB]	# o: flux error array
real	fwhm[ARB]	# o: wavelength error array
#--
int	irow, nrow
pointer	sp, units, unit2, wavval, fluxval, errval, fwhmval, nulflg
pointer	tp, cp
real	fillval

string	wavecol   "WAVELENGTH"
string	fluxcol   "FLUX"
string	errcol    "STATERROR"
string	fwhmcol   "FWHM"
string	nullwave  "Indef found in wavelength column"
string	badwave   "Wavelength units not found, angstroms assumed"
#string	badflux   "Flux units not found, PHOTLAM assumed"
string	badunits  "Units mismatch in spectrophotometry table"
string	nullthru  "Flux column all indef"
string	wavorder  "Wavelength column is not sorted"

bool	strne()
int	tbpsta(), anytoang()
int	fillnull(), is_magunit(), synsort4()
pointer	opnsyntab()

begin
	# Allocate memory for temporary string

	call smark (sp)
	call salloc (units, SZ_COLUNITS, TY_CHAR)
	call salloc (unit2, SZ_COLUNITS, TY_CHAR)

	# Open table containg spectrum

	tp = opnsyntab (table)

	# Allocate temporary arrays to hold contents of table

	nrow = tbpsta (tp, TBL_NROWS)
	call salloc (wavval, nrow, TY_REAL)
	call salloc (fluxval, nrow, TY_REAL)
	call salloc (errval, nrow, TY_REAL)
	call salloc (fwhmval, nrow, TY_REAL)
	call salloc (nulflg, nrow, TY_BOOL)

	# Read wavelength column

	call syncolptr (tp, wavecol, 1, cp)
	call tbcgtr (tp, cp, Memr[wavval], Memb[nulflg], 1, nrow)

	# Check for nulls

	do irow = 1, nrow {
	    if (Memb[nulflg+irow-1])
		call printerr_str (nullwave, table)
	}

	# Convert wavelength units to angstroms

	call tbcigt (cp, TBL_COL_UNITS, Memc[units], SZ_COLUNITS)

	if (anytoang (Memc[units], Memr[wavval], nrow) == NO) {
	    if (tbpsta (tp, TBL_WHTYPE) != TBL_TYPE_TEXT)
		call synphotwarn (badwave, table)
	}

	# Read fwhm column and remove nulls

	iferr {
	    call syncolptr (tp, fwhmcol, 4, cp)
	} then {
	    call amovkr (INDEFR, Memr[fwhmval], nrow)
	} else {
	    call tbcgtr (tp, cp, Memr[fwhmval], Memb[nulflg], 1, nrow)
	    call tbcigt (cp, TBL_COL_UNITS, Memc[unit2], SZ_COLUNITS)

	    if (strne (Memc[units], Memc[unit2]))
		call printerr_str (badunits, fwhmcol)

	    if (fillnull (0.0, nrow, Memr[wavval], Memr[fwhmval]) == ERR) {
		call amovkr (INDEFR, Memr[fwhmval], nrow)

	    } else {
		if (anytoang (Memc[units], Memr[fwhmval], nrow) == NO) {
		    if (tbpsta (tp, TBL_WHTYPE) != TBL_TYPE_TEXT)
			call synphotwarn (badwave, table)
		}
	    }
	}

	# Read flux column and remove nulls

	call syncolptr (tp, fluxcol, 2, cp)
	call tbcgtr (tp, cp, Memr[fluxval], Memb[nulflg], 1, nrow)
	call tbcigt (cp, TBL_COL_UNITS, Memc[units], SZ_COLUNITS)

	if (is_magunit (Memc[units]) == NO) {
	    fillval = 0.0
	} else {
	    fillval = 100.0
	}

	if (fillnull (fillval, nrow, Memr[wavval], Memr[fluxval]) == ERR)
	    call printerr_str (nullthru, table)

	# Read error column and remove nulls

	iferr {
	    call syncolptr (tp, errcol, 3, cp)
	} then {
	    call amovkr (INDEFR, Memr[errval], nrow)
	} else {
	    call tbcgtr (tp, cp, Memr[errval], Memb[nulflg], 1, nrow)
	    call tbcigt (cp, TBL_COL_UNITS, Memc[unit2], SZ_COLUNITS)

	    if (strne (Memc[units], Memc[unit2]))
		call printerr_str (badunits, errcol)

	    if (fillnull (0.0, nrow, Memr[wavval], Memr[errval]) == ERR)
		call amovkr (INDEFR, Memr[errval], nrow)
	}

	# Convert flux and error units to PHOTLAM

	call errphot (Memc[units], nrow, Memr[wavval], 
		      Memr[fluxval], Memr[errval])

	# Sort the arrays in ascending order and remove nulls

	# Sort arrays  on ascending order of wavelength
	# and interpolate on wavelength grid

	if (synsort4 (nrow, Memr[wavval], Memr[fluxval], 
		      Memr[errval], Memr[fwhmval]) == ERR)
	    call printerr_str (wavorder, table)

	call syninterp (nrow, Memr[wavval], Memr[fluxval], 
			nwave, wave, flux)

	if (IS_INDEFR(Memr[errval])) {
	    call amovkr (INDEFR, err, nwave)
	} else {
	    call syninterp (nrow, Memr[wavval], Memr[errval], 
			    nwave, wave, err)
	}

	if (IS_INDEFR(Memr[fwhmval])) {
	    call amovkr (INDEFR, fwhm, nwave)
	} else {
	    call syninterp (nrow, Memr[wavval], Memr[fwhmval], 
			    nwave, wave, fwhm)
	}

	call sfree (sp)
end

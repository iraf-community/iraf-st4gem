include	<tbset.h>

# RDBAND -- Read a bandpass from a table sampled on the wavelength set

procedure rdband (filename, nwave, wave, band)

char	filename[ARB]	# i: name of file containing bandpass
int     nwave           # i: length of wavelength and bandpass arrays
real    wave[ARB]       # i: wavelengths at which bandpass is computed
real	band[ARB]	# o: throughput at sampled wavelengths
#--
int	ic, nc, irow, nrow, thrunum
pointer	sp, thruname, tabname, units
pointer	tp, wv, thu, wavval, thruval, nulflg

string	wavecol   "WAVELENGTH"
string	thrucol   "THROUGHPUT"
string	nullwave  "Indef found in wavelength column"
string	badwave   "Wavelength units not found, angstroms assumed"
string	nullthru  "Throughput column all indef"
string	wavorder  "Wavelength column is not sorted"

int	ctoi(), tbpsta(), anytoang(), fillnull(), synsort2()
pointer	opnsyntab()

errchk	opnsyntab, syncolptr, synphoterr

begin
	# Allocate memory for temporary string

	call smark (sp)
	call salloc (thruname, SZ_FNAME, TY_CHAR)
	call salloc (tabname, SZ_FNAME, TY_CHAR)
	call salloc (units, SZ_COLUNITS, TY_CHAR)

	# Extract table and column name for file name and open table

	call breakcomp (filename, Memc[tabname], Memc[thruname], SZ_FNAME)
	tp = opnsyntab (Memc[tabname])

	# Allocate temporary arrays to hold contents of table

	nrow = tbpsta (tp, TBL_NROWS)
	call salloc (wavval, nrow, TY_REAL)
	call salloc (thruval, nrow, TY_REAL)
	call salloc (nulflg, nrow, TY_BOOL)

	# Get default throughput column name if not found in file name

	if (Memc[thruname] == EOS) {
	    thrunum = 2
	    call strcpy (thrucol, Memc[thruname], SZ_FNAME)

	} else {
	    ic = 1
	    nc = ctoi (Memc[thruname], ic, thrunum)
	}

	# Read throughput array and remove nulls

	call syncolptr (tp, Memc[thruname], thrunum, thu)
	call tbcgtr (tp, thu, Memr[thruval], Memb[nulflg], 1, nrow)

	if (fillnull (0.0, nrow, Memr[wavval], Memr[thruval]) == ERR)
	    call synphoterr (nullthru, filename)

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

	# Sort throughput and wavelength in ascending order 
	# and interpolate on wavelength grid

	if (synsort2 (nrow, Memr[wavval], Memr[thruval]) == ERR)
	    call synphoterr (wavorder, filename)

	call syninterp (nrow, Memr[wavval], Memr[thruval], 
			nwave, wave, band)

	call sfree (sp)
end

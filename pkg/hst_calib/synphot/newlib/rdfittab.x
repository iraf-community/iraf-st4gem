include	<tbset.h>

#* HISTORY *
#* B. Simon	29-Jul-94	original

procedure rdfittab (input, equal, indegree, fluxunits, 
		    ndata, wave, data, weight)

char	input[ARB]	# i: file containing observed data
bool	equal		# i: use equal weighting on the observed data?
int	indegree	# o: degree of the data in the input file
char	fluxunits[ARB]	# o: flux units of data column
int	ndata		# o: number of data pointes
pointer	wave		# o: wavelengths at which fit is evaluated
pointer	data		# o: observed data to be fit
pointer	weight		# o: weights used in computing residuals
#--
bool	same
int	nrow, done
pointer	tp, cp, sp, wavunits, colflag, tabflag

string	wavecol    "WAVELENGTH"
string	fluxcol    "FLUX"
string	fluxerrcol "STATERROR"
string	thrucol    "THROUGHPUT"
string	thruerrcol "ERROR"

string	notabcol   "THROUGHPUT or FLUX column not found in table"
string	nulldata   "No valid data in input file"
string	badwave    "Wavelength units not found, angstroms assumed"
string	wavorder   "Wavelength column is not sorted"

int	tbpsta(), anytoang(), anytophot(), synsort3()
pointer	tbtopn()

begin
	# Allocate temporary memory for strings

	call smark (sp)
	call salloc (wavunits, SZ_COLUNITS, TY_CHAR)

	# Open input table

	tp = tbtopn (input, READ_ONLY, NULL)
	nrow = tbpsta (tp, TBL_NROWS)

	# Allocate memory for output arrays

	call salloc (colflag, nrow, TY_BOOL)
	call salloc (tabflag, nrow, TY_BOOL)

	call malloc (wave, nrow, TY_REAL)
	call malloc (data, nrow, TY_REAL)
	call malloc (weight, nrow, TY_REAL)

	# Read wavelength column

	call syncolptr (tp, wavecol, 1, cp)
	call tbcgtr (tp, cp, Memr[wave], Memb[tabflag], 1, nrow)
	call tbcigt (cp, TBL_COL_UNITS, Memc[wavunits], SZ_COLUNITS)

	# Read data column. Use column name to determine table type

	iferr {
	    call syncolptr (tp, fluxcol, 2, cp)
	} then {
	    iferr {
		call syncolptr (tp, thrucol, 2, cp) 
	    } then {
		call printerr_str (notabcol, input)
	    } else {
		indegree = 0
		fluxunits[1] = EOS
		call tbcgtr (tp, cp, Memr[data], Memb[colflag], 1, nrow)
	    }
	} else {
	    indegree = 1
	    call tbcgtr (tp, cp, Memr[data], Memb[colflag], 1, nrow)
	    call tbcigt (cp, TBL_COL_UNITS, fluxunits, SZ_COLUNITS)
	}

	call alogor (Memb[tabflag], Memb[colflag], Memb[tabflag], nrow)

	# Read error column

	if (equal) {
	    same = true
	    call amovkr (1.0, Memr[weight], nrow)

	} else if (indegree == 0) {
	    iferr {
		call syncolptr (tp, thruerrcol, 3, cp)
	    } then {
		same = true
		call amovkr (1.0, Memr[weight], nrow)
	    } else {
		same = false
		call tbcgtr (tp, cp, Memr[weight], Memb[colflag], 1, nrow)
		call flagneg (nrow, Memr[weight], Memb[colflag])
		call alogor (Memb[tabflag], Memb[colflag], Memb[tabflag], nrow)
	    }

	} else {
	    iferr {
		call syncolptr (tp, fluxerrcol, 3, cp)
	    } then {
		same = true
		call amovkr (1.0, Memr[weight], nrow)
	    } else {
		same = false
		call tbcgtr (tp, cp, Memr[weight], Memb[colflag], 1, nrow)
		call flagneg (nrow, Memr[weight], Memb[colflag])
		call alogor (Memb[tabflag], Memb[colflag], Memb[tabflag], nrow)
	    }
	}

	# Remove nulls from columns

	call rmvnull (nrow, ndata, Memb[tabflag], Memr[wave])
	call rmvnull (nrow, ndata, Memb[tabflag], Memr[data])
	call rmvnull (nrow, ndata, Memb[tabflag], Memr[weight])

	if (ndata == 0)
	    call printerr_str (nulldata, input)

	# Convert wavelength units to angstroms

	if (anytoang (Memc[wavunits], Memr[wave], ndata) == NO) {
	    if (tbpsta (tp, TBL_WHTYPE) != TBL_TYPE_TEXT)
		call synphotwarn (badwave, input)
	}

	# Convert flux and its error to photlam

	if (indegree == 1) {
	    if (same) {
		done = anytophot (fluxunits, ndata, 
				  Memr[wave], Memr[data])
	    } else {
		call errphot (fluxunits, ndata, Memr[wave], 
			      Memr[data], Memr[weight])
	    }
	}

	# Convert errors to weights

	call cvtweight (ndata, Memr[weight])

	# Sort the data in order of increasing wavelengths

	if (synsort3 (ndata, Memr[wave], Memr[data], Memr[weight]) == ERR)
	    call printerr_str (wavorder, input)

	call tbtclo (tp)
	call sfree (sp)
end

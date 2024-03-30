include	<tbset.h>

#* HISTORY*
#* B.Simon	04-Mar-93	original
#* B.Simon	09-Apr-93	changed to call rdtabform
#* B.Simon	07-Jan-94	add input form
#* B.Simon	03-Jun-94	revised for new synphot library

# RDTABSPEC -- Read a spectrum from a table

procedure rdtabspec (table, inform, badpix, ilen, iwave, ispec)

char	table[ARB]	# i: table name
char	inform[ARB]	# i: spectral form string
real	badpix		# i: replacement value for bad pixel
int	ilen		# o: number of wavelengths
pointer	iwave		# o: wavelength array
pointer	ispec		# o: spectral array
#--
int	irow, done
pointer	tp, sp, form, wtype, wavptr, specptr, nulbuf

string	nowave    "Wavelength column not found in table"
string	nospec    "Flux column not found in table"
string	nullwave  "Indef found in wavelength column"

int	tbpsta(), fillnull(), anytoang(), anytophot()
pointer	tbtopn()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (form, SZ_FNAME, TY_CHAR)
	call salloc (wtype, SZ_FNAME, TY_CHAR)

	# Open table and get column pointers

	tp = tbtopn (table, READ_ONLY, 0)

	if (tbpsta (tp, TBL_WHTYPE) == TBL_TYPE_TEXT) {
	    call tbcfnd (tp, "c1", wavptr, 1)
	    call tbcfnd (tp, "c2", specptr, 1)

	} else {
	    call tbcfnd (tp, "WAVELENGTH", wavptr, 1)
	    call tbcfnd (tp, "FLUX", specptr, 1)
	    if (specptr == NULL)
		call tbcfnd (tp, "FLUX1", specptr, 1)
	}

	if (wavptr == NULL) {
	    call tbtclo (tp)
	    call printerr_str (nowave, table)
	}

	if (specptr == NULL) {
	    call tbtclo (tp)
	    call printerr_str (nospec, table)
	}

	# Read wavelength and flux columns

	ilen = tbpsta (tp, TBL_NROWS)
	call malloc (iwave, ilen, TY_REAL)
	call malloc (ispec, ilen, TY_REAL)
	call malloc (nulbuf, ilen, TY_BOOL)

	call tbcgtr (tp, wavptr, Memr[iwave], Memb[nulbuf], 1, ilen)
	call tbcgtr (tp, specptr, Memr[ispec], Memb[nulbuf], 1, ilen)

	# Read column units

	call tbcigt (wavptr, TBL_COL_UNITS, Memc[wtype], SZ_FNAME)
	call tbcigt (specptr, TBL_COL_UNITS, Memc[form], SZ_FNAME)

	if (Memc[form] == EOS) {
	    call strcpy (inform, Memc[form], SZ_FNAME)
	} else {
	    call strfix (Memc[form])
	}

	# Check for null wavelengths

	do irow = 1, ilen {
	    if (Memb[nulbuf+irow-1])
		call printerr_str (nullwave, table)
	}

	# Convert the units to angstroms and photlam

	done = fillnull (badpix, ilen, Memr[iwave], Memr[ispec])
	done = anytoang (Memc[wtype], Memr[iwave], ilen)
	done = anytophot (Memc[form], ilen, Memr[iwave], Memr[ispec])

	# Put spectrum in increasing monotonic order

	call synsort2 (ilen, Memr[iwave], Memr[ispec])

	call tbtclo (tp)
	call mfree (nulbuf, TY_BOOL)
	call sfree (sp)
end

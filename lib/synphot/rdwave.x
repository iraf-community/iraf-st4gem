include <tbset.h>

# RDWAVE -- Read wavelength array from table

procedure rdwave (tp, nwave, wave)

pointer	tp		# i: table descriptor
int	nwave		# i: number of wavelengths
real	wave[ARB]	# o: wavelength array
#--
int	iwave
pointer	wv, sp, nulflg, units, tabname

string	wavecol   "WAVELENGTH"
string	nullwave  "Indef found in wavelength table"
string	badtype   "Illegal data type for wavelength column"
string	badunits  "Wavelength units not found, angstroms assumed"
string	wavorder  "Wavelength table is not sorted"

int	tbpsta(), tbcigi(), anytoang(), synsort1()

errchk	tbcgtr, syncolptr, anytoang, synsort1

begin
	# Allocate memory for temporary array and strings

	call smark (sp)
	call salloc (nulflg, nwave, TY_BOOL)
	call salloc (units, SZ_COLUNITS, TY_CHAR)
	call salloc (tabname, SZ_FNAME, TY_CHAR)

	# Read wavelength array

	call syncolptr (tp, wavecol, 1, wv)
	call tbcgtr (tp, wv, wave, Memb[nulflg], 1, nwave)

	# Check for string column

	if (tbcigi (wv, TBL_COL_DATATYPE) < 0)
	    call synphoterr (badtype, "string")

	# Check for nulls

	do iwave = 1, nwave {
	    if (Memb[nulflg+iwave-1]) {
		call tbtnam (tp, Memc[tabname], SZ_FNAME)
		call synphoterr (nullwave, Memc[tabname])
	    }
	}

	# Convert units to angstroms

	call tbcigt (wv, TBL_COL_UNITS, Memc[units], SZ_COLUNITS)

	if (anytoang (Memc[units], wave, nwave) == NO) {
	    if (tbpsta (tp, TBL_WHTYPE) != TBL_TYPE_TEXT) {
		call tbtnam (tp, Memc[tabname], SZ_FNAME)
		call synphotwarn (badunits, Memc[tabname])
	    }
	}
	
	# Put in ascending order

	if (synsort1 (nwave, wave) == ERR) {
	    call tbtnam (tp, Memc[tabname], SZ_FNAME)
	    call synphoterr (wavorder, Memc[tabname])
	}

	call sfree (sp)
end

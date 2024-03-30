include <tbset.h>
include "libsynphot.h"

# EVALFILT -- Compute the througput of a component on a wavelength grid

procedure evalfilt (filename, nparam, param, nwave, wave, noerr, 
		    compthru, comperr)

char	filename[ARB]		# i: component filename
int	nparam			# i: number of component parameters
real	param[ARB]		# i: component parameters
int	nwave			# i: number of wavelengths
real	wave[ARB]		# i: wavelength grid
bool	noerr			# u: do not retrieve throughput error?
real	compthru[ARB]		# o: component throughput
real	comperr[ARB]		# o: component throughput error
#--
bool	reverse
char	pchar
int	nc, thrunum, errnum, irow, nrow
pointer	sp, thruname, errname, tabname, units
pointer	tp, wv, wavval, thruval, errval, nulflg

data	pchar	/ PCH /

string	wavecol   "WAVELENGTH"
string	thrucol   "THROUGHPUT"
string	errcol    "ERROR"
string	errsuffix "_err"
string	nullwave  "Indef found in wavelength column"
string	badunits  "Wavelength units not found, angstroms assumed"
string	nullthru  "Throughput column all indef"
string	wavorder  "Wavelength column is not sorted"

int	tbpsta(), stridx(), anytoang(), wavedir()
pointer	opnsyntab()

errchk	opnsyntab, syncolptr, synphoterr, anytoang
errchk	tbcigt, tbcgtr, synsort2, synsort3

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (thruname, SZ_FNAME, TY_CHAR)
	call salloc (errname, SZ_FNAME, TY_CHAR)
	call salloc (tabname, SZ_FNAME, TY_CHAR)
	call salloc (units, SZ_COLUNITS, TY_CHAR)

	# Extract table and column name for file name and open table

	call breakcomp (filename, Memc[tabname], Memc[thruname], SZ_FNAME)
	tp = opnsyntab (Memc[tabname])

	# Allocate memory for table columns

	nrow = tbpsta (tp, TBL_NROWS)
	call salloc (wavval, nrow, TY_REAL)
	call salloc (thruval, nrow, TY_REAL)
	call salloc (errval, nrow, TY_REAL)
	call salloc (nulflg, nrow, TY_BOOL)

	# Get default column nmaes if not found in file name
	# Otherwise, build error column name from throughput column name

	if (Memc[thruname] == EOS) {
	    thrunum = 2
	    errnum = 3
	    call strcpy (thrucol, Memc[thruname], SZ_FNAME)
	    call strcpy (errcol, Memc[errname], SZ_FNAME)

	} else {
	    thrunum = 0
	    errnum = 0

	    nc = stridx (pchar, Memc[thruname])
	    if (nc == 0) {
		call strcpy (Memc[thruname], Memc[errname], SZ_FNAME)
		call strcat (errsuffix, Memc[errname], SZ_FNAME)
	    } else {
		call strcpy (Memc[thruname], Memc[errname], nc)
		call strcat (errsuffix, Memc[errname], SZ_FNAME)
		call strcat (Memc[thruname+nc-1], Memc[errname], SZ_FNAME)
	    }
	}

	# Get wavelength column and check for nulls

	call syncolptr (tp, wavecol, 1, wv)
	call tbcigt (wv, TBL_COL_UNITS, Memc[units], SZ_COLUNITS)
	call tbcgtr (tp, wv, Memr[wavval], Memb[nulflg], 1, nrow)
	do irow = 1, nrow {
	    if (Memb[nulflg+irow-1])
		call synphoterr (nullwave, Memc[tabname])
	}

	# Convert wavelength units to angstroms

	if (anytoang (Memc[units], Memr[wavval], nrow) == NO) {
	    if (tbpsta (tp, TBL_WHTYPE) != TBL_TYPE_TEXT)
		call synphotwarn (badunits, Memc[tabname])
	}

	# Reverse wavelength order, if necessary

	switch (wavedir (nrow, Memr[wavval])) {
	case 1:
	    reverse = false
        case 0:
	    call synphoterr (wavorder, Memc[tabname])
	case -1:
	    reverse = true
	    call flip (nrow, Memr[wavval])
	}


	# Get throughput column

	call getthrucol (reverse, tp, Memc[thruname], thrunum, nparam, 
			 param, nrow, Memr[wavval], Memr[thruval])

	# Get throughput error. Error return means that 
	# column not found or all null

	if (! noerr) {
	    iferr {
		call getthrucol (reverse, tp, Memc[errname], errnum, nparam, 
				 param, nrow, Memr[wavval], Memr[errval])
	    } then {
		noerr = true
	    }
	}

	# Interpolate on wavelength grid. Separate cases depending 
	# on whether throughput error present or not

	if (noerr) {
	    call syninterp (nrow, Memr[wavval], Memr[thruval], 
			    nwave, wave, compthru)

	    call amovkr (INDEFR, comperr, nwave)

	} else {
	    call syninterp (nrow, Memr[wavval], Memr[thruval], 
			    nwave, wave, compthru)

	    call syninterp (nrow, Memr[wavval], Memr[errval], 
			    nwave, wave, comperr)
	}

	call sfree (sp)
end

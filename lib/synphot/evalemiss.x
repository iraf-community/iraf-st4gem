include <tbset.h>
include "libsynphot.h"

# EVALEMISS -- Compute the emission of a component on a wavelength grid
#   This routine was cloned from evalfilt.x - V. G. Laidler, sep 2001

procedure evalemiss (filename, nparam, param, nwave, wave, noerr, 
		    compemiss, comperr)

char	filename[ARB]		# i: component filename
int	nparam			# i: number of component parameters
real	param[ARB]		# i: component parameters
int	nwave			# i: number of wavelengths
real	wave[ARB]		# i: wavelength grid
bool	noerr			# u: do not retrieve throughput error?
real	compemiss[ARB]		# o: component throughput
real	comperr[ARB]		# o: component throughput error
#--
bool	reverse
char	pchar
int	nc, emissnum, errnum, irow, nrow
pointer	sp, emissname, errname, tabname, units
pointer	tp, wv, wavval, emissval, errval, nulflg

data	pchar	/ PCH /

string	wavecol   "WAVELENGTH"
string	colname   "EMISSIVITY"
string	errcol    "ERROR"
string	errsuffix "_err"
string	nullwave  "Indef found in wavelength column"
string	badunits  "Wavelength units not found, angstroms assumed"
string	wavorder  "Wavelength column is not sorted"

int	tbpsta(), stridx(), anytoang(), wavedir()
pointer	opnsyntab()

errchk	opnsyntab, syncolptr, synphoterr, anytoang
errchk	tbcigt, tbcgtr, synsort2, synsort3

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (emissname, SZ_FNAME, TY_CHAR)
	call salloc (errname, SZ_FNAME, TY_CHAR)
	call salloc (tabname, SZ_FNAME, TY_CHAR)
	call salloc (units, SZ_COLUNITS, TY_CHAR)

	# Extract table and column name for file name and open table

	call breakcomp (filename, Memc[tabname], Memc[emissname], SZ_FNAME)
	tp = opnsyntab (Memc[tabname])

	# Allocate memory for table columns

	nrow = tbpsta (tp, TBL_NROWS)
	call salloc (wavval, nrow, TY_REAL)
	call salloc (emissval, nrow, TY_REAL)
	call salloc (errval, nrow, TY_REAL)
	call salloc (nulflg, nrow, TY_BOOL)

	# Get default column nmaes if not found in file name
	# Otherwise, build error column name from throughput column name

	if (Memc[emissname] == EOS) {
	    emissnum = 2
	    errnum = 3
	    call strcpy (colname, Memc[emissname], SZ_FNAME)
	    call strcpy (errcol, Memc[errname], SZ_FNAME)

	} else {
	    emissnum = 0
	    errnum = 0

	    nc = stridx (pchar, Memc[emissname])
	    if (nc == 0) {
		call strcpy (Memc[emissname], Memc[errname], SZ_FNAME)
		call strcat (errsuffix, Memc[errname], SZ_FNAME)
	    } else {
		call strcpy (Memc[emissname], Memc[errname], nc)
		call strcat (errsuffix, Memc[errname], SZ_FNAME)
		call strcat (Memc[emissname+nc-1], Memc[errname], SZ_FNAME)
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


	# Get emissivity column.
	# Note that in spite of the subroutine name, getthrucol is a perfectly general
	# routine (except for the text of its error messages), and can be used to read emissivity 
        # data given the proper calling sequence. VGL

	call getthrucol (reverse, tp, Memc[emissname], emissnum, nparam, 
			 param, nrow, Memr[wavval], Memr[emissval])

	# Get emissivity error. Error return means that 
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
	    call syninterp (nrow, Memr[wavval], Memr[emissval], 
			    nwave, wave, compemiss)

	    call amovkr (INDEFR, comperr, nwave)

	} else {
	    call syninterp (nrow, Memr[wavval], Memr[emissval], 
			    nwave, wave, compemiss)

	    call syninterp (nrow, Memr[wavval], Memr[errval], 
			    nwave, wave, comperr)
	}


	call sfree (sp)
end

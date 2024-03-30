include	<tbset.h>
include "libsynphot.h"

#* HISTORY*
#* B.Simon	31-Oct-94	Use different null values based on col units
#* B.Simon	20-Oct-97	Allow for no thriughput column

# WAVELIMITS -- Compute the extreme wavelengths for a throughput column

procedure wavelimits (tp, wv, thu, wlow, whi)

pointer	tp		# i: table descriptor
pointer	wv		# i: wavelength column
pointer	thu		# i: throughput column pointer
real	wlow		# o: low wavelength
real	whi		# o: high wavelength
#--
int	irow, jrow, nrow
pointer	sp, tabname, units
real	nullval, trans, oldtrans, w[2]

string	nullwave  "Null wavelength in throughput table"
string	badunits  "Wavelength units not found, angstroms assumed"

int	tbpsta(), is_magunit(), anytoang()
errchk	synphoterr, anytoang, tbegtr

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (tabname, SZ_FNAME, TY_CHAR)
	call salloc (units, SZ_COLUNITS, TY_CHAR)

	# Set the null value according to the throughput column type

	if (thu != NO_COLUMN) {
	    call tbcigt (thu, TBL_COL_UNITS, Memc[units], SZ_COLUNITS)
	    if (is_magunit (Memc[units]) == NO) {
		nullval = 0.0
	    } else {
		nullval = 100.0
	    }
	}

	# Read number of rows and wavelength column units

	nrow = tbpsta (tp, TBL_NROWS)
	call tbcigt (wv, TBL_COL_UNITS, Memc[units], SZ_COLUNITS)

	# Find first wavelength in table associated with significant throughput

	jrow = 0
	call tbegtr (tp, wv, 1, w[1])
	if (thu != NO_COLUMN) {
	    call tbegtr (tp, thu, 1, oldtrans)

	    do irow = 2, nrow {
		call tbegtr (tp, thu, irow, trans)
		if (! IS_INDEFR(trans)) {
		    if (trans != nullval) {
			if (! IS_INDEFR(oldtrans)) {
			    jrow = irow - 1
			} else {
			    jrow = irow
			}
			call tbegtr (tp, wv, jrow, w[1])
			break
		    }
		}
		oldtrans = trans
	    }
	}

	# Find last wavelength in table associated with significant throughput

	jrow = 0
	call tbegtr (tp, wv, nrow, w[2])
	if (thu != NO_COLUMN) {
	    call tbegtr (tp, thu, nrow, oldtrans)

	    do irow = nrow-1, 1, -1 {
		call tbegtr (tp, thu, irow, trans)
		if (! IS_INDEFR(trans)) {
		    if (trans != nullval) {
			if (! IS_INDEFR(oldtrans)) {
			    jrow = irow + 1
			} else {
			    jrow = irow
			}
			call tbegtr (tp, wv, jrow, w[2])
			break
		    }
		}
		oldtrans = trans
	    }
	}

	# Check for indefinite wavelengths, which are not allowed

	if (IS_INDEFR(w[1]) || IS_INDEFR(w[2])) {
	    call tbtnam (tp, Memc[tabname], SZ_FNAME)
	    call synphoterr (nullwave, Memc[tabname])
	}

	# Convert wavelength units to angstroms

	if (anytoang (Memc[units], w, 2) == NO) {
	    if (tbpsta (tp, TBL_WHTYPE) != TBL_TYPE_TEXT) {
		call tbtnam (tp, Memc[tabname], SZ_FNAME)
		call synphotwarn (badunits, Memc[tabname])
	    }
	}

	# Copy wavelengths into outputs

	if (w[1] < w[2]) {
	    wlow = w[1]
	    whi = w[2]
	} else {
	    wlow = w[2]
	    whi = w[1]
	}

	call sfree (sp)
end

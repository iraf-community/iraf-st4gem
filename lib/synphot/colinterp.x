include <tbset.h>

# COL_INTERP -- Interpolate between columns in a parameterized table

procedure col_interp (reverse, tp, nparam, param, collist, nrow, wave, value)

bool	reverse		# i: reverse columns before interpolating?
pointer	tp		# i: table descriptor
int	nparam		# i: number of parameters
real	param[ARB]	# i: parameter values
pointer	collist[ARB]	# u: list of columns to interpolate between
int	nrow		# i: number of rows to read
real	wave[ARB]	# i: wavelength array
real	value[ARB]	# o: column values
#--
bool	extrap, shift
int	ic, nc, icol, jcol, ncol, ipar, nval, status
pointer	sp, nulflg, parval, partype, params, tabname, parname
pointer	checked, thruval, minthru, maxthru, newthru
real	loval, hival

string	outside   "Extrapolating in parameterized throughput table"
string	nulldata  "Table column values are all INDEF"

bool	streq()
int	fillnull(), word_fetch()

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (nulflg, nrow, TY_BOOL)
	call salloc (parval, nparam, TY_REAL)
	call salloc (partype, SZ_FNAME, TY_CHAR)
	call salloc (params, SZ_FNAME, TY_CHAR)
	call salloc (tabname, SZ_FNAME, TY_CHAR)
	call salloc (parname, SZ_COLNAME, TY_CHAR)

	ncol = 2 ** nparam
	call salloc (checked, ncol, TY_INT)
	call salloc (thruval, ncol, TY_POINTER)

	# Check for null values in column list and replace them 
	# with their paired values

	extrap = false
	for (jcol = 1; jcol < ncol; jcol = 2 * jcol) {
	    call amovki (NO, Memi[checked], ncol)

	    do icol = 1, ncol-jcol {
		if (Memi[checked+icol-1] == NO) {
		    if (collist[icol] == NULL) {
			collist[icol] = collist[icol+jcol]
			extrap = true
		    }

		    if (collist[icol+jcol] == NULL) {
			collist[icol+jcol] = collist[icol]
			extrap = true
		    }

		    Memi[checked+icol+jcol-1] = YES
		}
	    }
	}

	if (extrap) {
	    call tbtnam (tp, Memc[tabname], SZ_FNAME)
	    call synphotwarn (outside, Memc[tabname])
	}

	# Read columns that bracket the value we are interpolating for

	do icol = 1, ncol {
	    call salloc (Memi[thruval+icol-1], nrow, TY_REAL)

	    call tbcgtr (tp, collist[icol], Memr[Memi[thruval+icol-1]],
			 Memb[nulflg], 1, nrow)

	    status = fillnull (0.0, nrow, wave, Memr[Memi[thruval+icol-1]])

	    if (status == ERR) {
		call tbtnam (tp, Memc[tabname], SZ_FNAME)
		call synphoterr (nulldata, Memc[tabname])
	    }

	    if (reverse)
		call flip (nrow, Memr[Memi[thruval+icol-1]])
	}

	# Read the header keyword that indicates if the interpolation
	# includes a shift of origin

	iferr {
	    call rdtabhdt (tp, "PARAMS", Memc[params], SZ_FNAME)
	} then {
	    Memc[params] = EOS
	}

	# Interpolate between successive arrays of column data

	ic = 1
	do ipar = 1, nparam {
	    icol = 1
	    jcol = 1

	    nc = word_fetch (Memc[params], ic, Memc[partype], SZ_FNAME)
	    call strlwr (Memc[partype])

	    shift = streq (Memc[partype], "wavelength")

	    while (icol < ncol) {
		call tbcigt (collist[icol], TBL_COL_NAME, 
			     Memc[parname], SZ_COLNAME)
		call breakparam (Memc[parname], nparam, nval, Memr[parval])
		loval = Memr[parval+ipar-1]

		call tbcigt (collist[icol+1], TBL_COL_NAME, 
			     Memc[parname], SZ_COLNAME)
		call breakparam (Memc[parname], nparam, nval, Memr[parval])
		hival = Memr[parval+ipar-1]

		minthru = Memi[thruval+icol-1]
		maxthru = Memi[thruval+icol]

		newthru = Memi[thruval+jcol-1]
		collist[jcol] = collist[icol]

		call shiftinterp (shift, loval, hival, param[ipar], 
				  nrow, wave, Memr[minthru], 
				  Memr[maxthru], Memr[newthru])

		icol = icol + 2
		jcol = jcol + 1
	    }

	    ncol = jcol - 1
	}

	# Copy interpolated values to output array

	call amovr (Memr[Memi[thruval]], value, nrow)
	call sfree (sp)
end

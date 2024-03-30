include <tbset.h>

# GETTHRUCOL -- Get the throughput values from throughput column(s)

procedure getthrucol (reverse, tp, colname, colnum, nparam, param, 
		      nrow, wave, value)

bool	reverse		# i: reverse value column?
pointer	tp		# i: table descriptor
char	colname[ARB]	# i: table column name
int	colnum		# i: column number
int	nparam		# i: number of parameters
real	param[ARB]	# i: parameter values
int	nrow		# i: number of rows to read
real	wave[ARB]	# u: wavelength array
real	value[ARB]	# o: column values
#--
int	ncol, nptr, iptr, utype, status
pointer	sp, nulflg, tabname, units, colptr

string	thrudic   "|transmission|qe|dn/photon|emissivity"
string	badunits  "Throughput column has illegal units"
string	noparams  "No parameters were given for parameterized component"
string	nulldata  "Table column values are all INDEF"

int	tbpsta(), strdic(), fillnull()

errchk	tbcgtr, thruputcol, getparcol, synphoterr

begin
	# Allocate meory for strings and temporary arrays

	call smark (sp)
	call salloc (nulflg, nrow, TY_BOOL)
	call salloc (tabname, SZ_FNAME, TY_CHAR)
	call salloc (units, SZ_COLUNITS, TY_CHAR)

	ncol = tbpsta (tp, TBL_NCOLS)
	call salloc (colptr, ncol, TY_POINTER)

	# Get columns that match column name (more than one if parameterized)

	call thruputcol (tp, colname, colnum, ncol, nptr, Memi[colptr])

	# Check throughput column units

	do iptr = 1, nptr {
	    call tbcigt (Memi[colptr], TBL_COL_UNITS, Memc[units], SZ_COLUNITS)

	    call strfix (Memc[units])
	    if (Memc[units] != EOS) {
		utype = strdic (Memc[units], Memc[units], SZ_COLUNITS, thrudic)
		if (utype == 0) {
		    call tbtnam (tp, Memc[tabname], SZ_FNAME)
		    call synphoterr (badunits, Memc[tabname])
		}
	    }
	}

	# Read column if not parameterized, otherwise interpolate between
	# parameterized columns that bracket parameters

	if (nparam == 0) {
	    if (nptr > 1)
		call synphoterr (noparams, colname)

	    call tbcgtr (tp, Memi[colptr], value, Memb[nulflg], 1, nrow)

	    status = fillnull (0.0, nrow, wave, value)
	    if (status == ERR) {
		call tbtnam (tp, Memc[tabname], SZ_FNAME)
		call synphoterr (nulldata, Memc[tabname])
	    }

	    if (reverse)
		call flip (nrow, value)

	} else {
	    call getparcol (reverse, tp, nptr, Memi[colptr], colname, 
			    nparam, param, nrow, wave, value)
	}

	call sfree (sp)
end

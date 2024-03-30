include	<tbset.h>
include "libsynphot.h"

# TABRANGE -- Compute the wavelength range of a tabulated passband or spectrum

procedure tabrange (fname, minwave, maxwave, status)

char	fname[ARB]	# i: table name
real	minwave		# o: low wavelength
real	maxwave		# o: high wavelength
int	status		# i: OK if wavelengths could be read
#--
char	pchar
pointer	sp, tabname, thruname, tp, wv, dep

data	pchar	 / PCH /	# character indicating parameterized column

string	wavecol  "WAVELENGTH"
string	fluxcol  "FLUX"
string	thrucol  "THROUGHPUT"

int	tbcigi(), stridx()
pointer	opnsyntab()
errchk	syncolptr, wavelimits

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (tabname, SZ_FNAME, TY_CHAR)
	call salloc (thruname, SZ_FNAME, TY_CHAR)

	# Read wavelength column

	call breakcomp (fname, Memc[tabname], Memc[thruname], SZ_FNAME)

	tp = opnsyntab (Memc[tabname])
	call syncolptr (tp, wavecol, 1, wv)

	# If column is type string, this is a filename list
	# set wavelength range to INDEF and set an error condition
	# Otherwise, read the dependent column and get the limits

	if (tbcigi (wv, TBL_COL_DATATYPE) < 0) {
	    minwave = INDEFR
	    maxwave = INDEFR
	    status = ERR

	} else {
	    # Distinguish between a spectrum and throughput table by the
	    # presence of the throughput column

	    if (Memc[thruname] == EOS)
		call strcpy (thrucol, Memc[thruname], SZ_FNAME)

	    # Check for parameterized column name. No dependent
	    # column is used when getting the range limits for
	    # a parameterized column, since we do not know
	    # which columns in the table will be used

	    if (stridx (pchar, Memc[thruname]) == 0) {
		iferr {
		    call syncolptr (tp, Memc[thruname], 0, dep)
		} then {
		    call syncolptr (tp, fluxcol, 2, dep)
		}

	    } else {
		dep = NO_COLUMN
	    }

	    call wavelimits (tp, wv, dep, minwave, maxwave)
	    status = OK
	}

	call sfree (sp)

end

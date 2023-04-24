include	<time.h>
include	"colnames.h"

# OPENDB -- Opens output database table.
#
# Modified 2/22/96 (I Busko) to implement chi-squared minimization.
# Modified 9/14/98 (IB) to write dummy row at table creation time. This
#    enables later creation of new columns in FITS tables.

pointer procedure opendb (outname)

char	outname[ARB]			# output table name

pointer	c, c1, tp
int	iomode

int	strlen(), tbtopn(), tbtacc()

errchk	tbtopen, tbtclo, tbcdef

begin
	if (strlen (outname) > 0) {
	    if (tbtacc(outname) == NO)
	        iomode = NEW_FILE
	    else
	        iomode = READ_WRITE

	    tp = tbtopn (outname, iomode, 0)

	    # create new table
	    if (iomode == NEW_FILE) {
	        call tbcdef (tp, c1, DB_CFILE, DB_UFILE, DB_FFILE, -DB_SFILE, 1, 1)
	        call tbcdef (tp, c, DB_CTIME, DB_UTIME, DB_FTIME, -DB_STIME, 1, 1)
	        call tbcdef (tp, c, DB_CFUNC, DB_UFUNC, DB_FFUNC, -DB_SFUNC, 1, 1)
	        call tbcdef (tp, c, DB_CUNIT, DB_UUNIT, DB_FUNIT, -DB_SUNIT, 1, 1)
	        call tbcdef (tp, c, DB_CDEGR, DB_UDEGR, DB_FDEGR, TY_INT, 1, 1)
	        call tbcdef (tp, c, DB_CNPTS, DB_UNPTS, DB_FNPTS, TY_INT, 1, 1)
	        call tbcdef (tp, c, DB_CCHI,  DB_UCHI,  DB_FCHI,  TY_REAL, 1, 1)
	        call tbcdef (tp, c, DB_CRMS,  DB_URMS,  DB_FRMS,  TY_REAL, 1, 1)
	        call tbcdef (tp, c, DB_CXMIN, DB_UXMIN, DB_FXMIN, TY_REAL, 1, 1)
	        call tbcdef (tp, c, DB_CXMAX, DB_UXMAX, DB_FXMAX, TY_REAL, 1, 1)
	        call tbtcre (tp)
		call tbeptt (tp, c1, 1, DUMMY)

	    }
	} else {
	    call error (0, "No output specification.")
	}

	return (tp)
end

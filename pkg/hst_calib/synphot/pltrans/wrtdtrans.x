include	<tbset.h>
include "dtrans.h"

#* HISTORY *
#* B.Simon	22-Jul-94	original

# WRTDTRANS -- Write results in a form readable by pltrans

procedure wrtdtrans (output, xmode, ymode, xform, yform, append, dtrans)

char	output[ARB]	# i: output table containing photometric data
char	xmode[ARB]	# i: x axis mode
char	ymode[ARB]	# i: y axis mode
char	xform[ARB]	# i: x axis form
char	yform[ARB]	# i: y axis form
bool	append		# i: append results to old table?
pointer	dtrans		# i: transformation data descriptor
#--
int	irow, neff
pointer	tp, cx, cy, xeff, yeff

string	blank  ""

int	tbtacc(), tbpsta()
pointer	tbtopn()

begin
	# Early exit if not output table or nothing to write

	if (output[1] == EOS || dtrans == NULL)
	    return

	# Open existing table if append mode
	# Otherwise create new table
 
	if (append && tbtacc (output) == YES) {
	    tp = tbtopn (output, READ_WRITE, NULL)
	    call syncolptr (tp, "FLUX1", 1, cx)
	    call syncolptr (tp, "FLUX2", 2, cy)

	    irow = tbpsta (tp, TBL_NROWS) + 1

	} else {
	    tp = tbtopn (output, NEW_FILE, NULL)
	    call tbpset (tp, TBL_MAXCOLS, 2)

	    call tbcdef (tp, cx, "FLUX1", xform, blank, TY_REAL, 1, 1)
	    call tbcdef (tp, cy, "FLUX2", yform, blank, TY_REAL, 1, 1)
	    call tbtcre (tp)

	    call tbhadt (tp, "XMODE", xmode)
	    call tbhadt (tp, "XFORM", xform)
	    call tbhadt (tp, "YMODE", ymode)
	    call tbhadt (tp, "YFORM", yform)

	    irow = 1
	}

	# Wtire data from structure to table columns

	neff = TRN_NEFF(dtrans)
	xeff = TRN_XEFF(dtrans)
	yeff = TRN_YEFF(dtrans)

	call tbcptr (tp, cx, Memr[xeff], irow, irow+neff-1)
	call tbcptr (tp, cy, Memr[yeff], irow, irow+neff-1)
	call tbtclo (tp)
end

include	<tbset.h>
include	"dtrans.h"

#* HISTORY *
#* B.Simon	21-Jul-94	original

# INPDTRANS -- Read photometric data from a previous run of pltrans

procedure inpdtrans (input, xmode, ymode, xform, yform, maxch, dtrans)

char	input[ARB]	# i: input table containing photometric data
char	xmode[ARB]	# u: x axis mode
char	ymode[ARB]	# u: y axis mode
char	xform[ARB]	# u: x axis form
char	yform[ARB]	# u: y axis form
int	maxch		# i: length of strings
pointer	dtrans		# o: transformation data descriptor
#--
int	nrow, irow, jrow
pointer	sp, tp, xflag, yflag, xeff, yeff, cx, cy

int	tbpsta()
pointer	tbtopn()

begin
	# Take early exit if no input table

	dtrans = NULL
	if (input[1] == EOS)
	    return

	# Open table and get column pointers

	tp = tbtopn (input, READ_ONLY, NULL)
	nrow = tbpsta (tp, TBL_NROWS)

	call syncolptr (tp, "FLUX1", 1, cx)
	call syncolptr (tp, "FLUX2", 2, cy)

	# Allocate arrays to hold table data

	call smark (sp)
	call salloc (xflag, nrow, TY_BOOL) 
	call salloc (yflag, nrow, TY_BOOL) 

	call malloc (dtrans, LEN_TRNSTRUCT, TY_STRUCT)
	call malloc (xeff, nrow, TY_REAL)
	call malloc (yeff, nrow, TY_REAL)

	TRN_NEFF(dtrans) = 0
	TRN_XEFF(dtrans) = xeff
	TRN_YEFF(dtrans) = yeff

	# Read header parameters and column units to
	# get values for mode and form that are possibly
	# different from those input as task parameters

	call hpovwrite (tp, "XMODE", xmode, maxch)
	call hpovwrite (tp, "XFORM", xform, maxch)
	call hpovwrite (tp, "YMODE", ymode, maxch)
	call hpovwrite (tp, "YFORM", yform, maxch)

	# Read photometric data from table

	call tbcgtr (tp, cx, Memr[xeff], Memb[xflag], 1, nrow)
	call tbcgtr (tp, cy, Memr[yeff], Memb[yflag], 1, nrow)

	# Remove null data from array and set length of arrays

	jrow = 1
	do irow = 1, nrow {
	    if (! Memb[xflag] && ! Memb[yflag]) {
		if (jrow < irow) {
		    Memr[xeff+jrow-1] = Memr[xeff+irow-1]
		    Memr[yeff+jrow-1] = Memr[yeff+irow-1]
		}

		jrow = jrow + 1
	    }
	}

	TRN_NEFF(dtrans) = jrow - 1

	# Clean up

	call tbtclo (tp) 
	call sfree (sp) 

end

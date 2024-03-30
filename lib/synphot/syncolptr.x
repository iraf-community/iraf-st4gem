include	<tbset.h>
include	"libsynphot.h"

# SYNCOLPTR -- Get a column pointer to a synphot table

procedure syncolptr (tp, colname, colnum, ptr)

pointer	tp		# i: table descriptor
char	colname[ARB]	# i: column name
int	colnum		# i: column number
pointer	ptr		# o: table column pointer
#--
string	nocolumn  "Column not found in table"

int	tbpsta()
pointer	tbcnum()
errchk synphoterr

begin
	# Retrieve columns for text tables by column number,
	# from binary tables by column name. To prevent the
	# use of text tables, set the column number to zero.

	if (tbpsta (tp, TBL_WHTYPE) == TBL_TYPE_TEXT) {
	    iferr (ptr = tbcnum (tp, colnum))
		ptr = NULL

	} else {
	    iferr (call tbcfnd (tp, colname, ptr, 1))
		ptr = NULL
	}

	if (ptr == NULL)
	    call synerrflag (COL_ERR, nocolumn, colname)

end

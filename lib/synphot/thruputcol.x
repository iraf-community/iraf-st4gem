include	<tbset.h>
include "libsynphot.h"

# THRUPUTCOL -- Return a list of throughput column pointers

procedure thruputcol (tp, colname, colnum, ncol, nptr, colptr)

pointer	tp		# i: table descriptor
char	colname[ARB]	# i: column name
int	colnum		# i: column number, for text table
int	ncol		# i: maximum number of table columns
int	nptr		# o: actual number of columns
pointer	colptr[ARB]	# o: array of column pointers
#--
char	pchar
int	nc, ic, icol, oldpar, newpar, newcol
pointer	sp, cp, oldname, newname

data	pchar	/ PCH /

string	nocolumn  "Column not found in table"
string	overcol   "Too many throughputs for parameterized component"

int	stridx(), tbpsta(), strncmp()
pointer	tbcnum()

errchk	syncolptr, synphoterr

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (oldname, SZ_COLNAME, TY_CHAR)
	call salloc (newname, SZ_COLNAME, TY_CHAR)

	nc = stridx (pchar, colname)

	if (nc == 0) {
	    # If the column name does not contain the parameter charater, 
	    # get the column pointer directly with syncolptr. 

	    call syncolptr (tp, colname, colnum, colptr[1])
	    nptr = 1

	} else {
	    # Convert column name to lower case for case insensitive match

	    call strcpy (colname, Memc[oldname], SZ_COLNAME)
	    call strlwr (Memc[oldname])

	    # Count number of parameters in input column name

	    oldpar = 0
	    for (ic = nc; colname[ic] != EOS; ic = ic + 1) {
		if (colname[ic] == pchar)
		    oldpar = oldpar + 1
	    }

	    # Compare each column name in the table with the input column 
	    # name to see if the leading portion and number of parameters 
	    # agree. If they do, copy the column pointer to the output array.

	    nptr = 0
	    newcol = tbpsta (tp, TBL_NCOLS)

	    do icol = 1, newcol {
		cp = tbcnum (tp, icol)

		call tbcigt (cp, TBL_COL_NAME, Memc[newname], SZ_COLNAME)
		call strlwr (Memc[newname])

		if (strncmp (Memc[oldname], Memc[newname], nc) == 0) {
		    newpar = 0
		    for (ic = nc; Memc[newname+ic-1] != EOS; ic = ic + 1) {
			if (Memc[newname+ic-1] == pchar)
			    newpar = newpar + 1
		    }

		    if (oldpar == newpar) {
			nptr = nptr + 1
			if (nptr > ncol)
			    call synphoterr (overcol, colname)

			colptr[nptr] = cp
		    }
		}
	    }

	    if (nptr == 0)
		call synphoterr (nocolumn, colname)
	}

	call sfree (sp)
end

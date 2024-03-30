include <tbset.h>
include <ctype.h>

# OPEN_COPY -- Open output table as possibly modified copy of input

pointer procedure open_copy (ofile, itp, title)

char	ofile[ARB]	# i: output file name
pointer	itp		# i: input table pointer
bool	title		# i: use first row as title?
#--
int	ic, nc, icol, ncol, colnum, datatype, lendata, lenfmt, index
pointer	otp, sp, colptr, colname, colunits, colfmt

int	tbpsta(), ctoi(), strlen(), word_find()
pointer	tbtopn(), tbcnum()

string	collist   "WAVELENGTH,THROUGHPUT,ERROR"
string	unitlist  "ANGSTROMS,TRANSMISSION,TRANSMISSION"
string	fmtlist   "%10.1f,%12.5g,%12.5g"

begin
	# Allocate temporary strings

	call smark (sp)
	call salloc (colname, SZ_COLNAME, TY_CHAR)
	call salloc (colunits, SZ_COLUNITS, TY_CHAR)
	call salloc (colfmt, SZ_COLFMT, TY_CHAR)

	otp = tbtopn (ofile, NEW_FILE, NULL)

	ncol = tbpsta (itp, TBL_NCOLS)
	do icol = 1, ncol {

	    # Get column information from input table

	    colptr = tbcnum (itp, icol)
	    call tbcinf (colptr, colnum, Memc[colname], Memc[colunits], 
			 Memc[colfmt], datatype, lendata, lenfmt)

	    # Set column names and units if ascii table

	    if (Memc[colname] == 'c' && IS_DIGIT (Memc[colname+1])) {
		ic = 1
		nc = strlen (Memc[colname+1])

		if (ctoi (Memc[colname+1], ic, index) == nc) {
		    if (title) {
			# Get column names from first row

			call tbegtt (itp, colptr, 1, Memc[colname],
				     SZ_COLNAME)

			nc = word_find (index, unitlist, Memc[colunits], 
					SZ_COLUNITS)

			nc = word_find (index, fmtlist, Memc[colfmt], 
					SZ_COLFMT)

		    } else if (ncol <= 3) {
			# Get column names from default values

			nc = word_find (index, collist, Memc[colname], 
					SZ_COLNAME)

			nc = word_find (index, unitlist, Memc[colunits], 
					SZ_COLUNITS)

			nc = word_find (index, fmtlist, Memc[colfmt], 
					SZ_COLFMT)
		    }
		}
	    }

	    # Copy column information to output table

	    datatype = TY_REAL
	    call tbcdef (otp, colptr, Memc[colname], Memc[colunits], 
			 Memc[colfmt], datatype, lendata, 1)
	}

	call tbtcre (otp)
	call sfree(sp)
	return (otp)
end


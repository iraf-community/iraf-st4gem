include	<tbset.h>

# COL_CATEGORY -- Assign a parameterized column to a catagory

procedure col_category (cp, nparam, param, ncat, category)

pointer	cp		# i: Column pointer
int	nparam		# i: Dimensionality of interpolating value
real	param[ARB]	# i: Interpolating value
int	ncat		# o: Number of categories column falls in
int	category[ARB]	# o: Categories column falls in
#--
int	nval, catnum, ipar, icat
pointer	sp, parval, parname

string	badname   "Number of parameters does not match column name"

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (parval, nparam, TY_REAL)
	call salloc (parname, SZ_COLNAME, TY_CHAR)

	# Get column name and extract parameters

	call tbcigt (cp, TBL_COL_NAME, Memc[parname], SZ_COLNAME)
	call breakparam (Memc[parname], nparam, nval, Memr[parval])

	if (nval != nparam)
	    call synphoterr (badname, Memc[parname]) 

	# Compare column parameters to interpolating parameters and
	# assign column to categories on basis of comparison.
	# Categories indicate whether the column parameter value is above 
	# or below the value we are interpolating for.  If the column 
	# parameter is equal to the value interpolated for, it is placed in 
	# both categories. 

	ncat = 1
	catnum = 1
	category[1] = 1

	# If the colum parameter is less than the interpolating parameter 
	# zero (implicitly) gets added to all values in the list.
	# If it is greater catnum gets added. If it is equal the list
	# is doubled in length and zero gets added to half the values
	# and catnum to the other half.

	do ipar = 1, nparam {
	    if (Memr[parval+ipar-1] > param[ipar]) {
		do icat = 1, ncat
		    category[icat] = category[icat] + catnum

	    } else if (Memr[parval+ipar-1] == param[ipar]) {
		do icat = ncat, 1, -1 {
		    category[2*icat] = category[icat] + catnum
		    category[2*icat-1] = category[icat]
		}

		ncat = 2 * ncat
	    }

	    catnum = 2 * catnum
	}

	call sfree (sp)
end

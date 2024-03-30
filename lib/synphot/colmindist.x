include <tbset.h>

# COL_MINDIST -- Find the columns in each category that bracket the parameter

procedure col_mindist (cp, nparam, param, ncat, category, 
		       scaledist, mindist, collist)

pointer	cp		# i: Column pointer
int	nparam		# i: Dimensionality of interpolating value
real	param[ARB]	# i: Interpolating value
int	ncat		# i: Number of categories column falls in
int	category[ARB]	# i: Categories column falls in
real	scaledist[ARB]	# i: Scaling factor used in distance calculation
real	mindist[ARB]	# u: Minimum distance so far
pointer	collist[ARB]	# u: List of column pointers at minimum distance
#--
int	nval, ipar, icat, jcat
pointer	sp, parval, parname
real	dist

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (parval, nparam, TY_REAL)
	call salloc (parname, SZ_COLNAME, TY_CHAR)

	# Get column name and extract parameters

	call tbcigt (cp, TBL_COL_NAME, Memc[parname], SZ_COLNAME)
	call breakparam (Memc[parname], nparam, nval, Memr[parval])

	# Compute distance between column parameter and interpolating value

	dist = 0.0
	do ipar = 1, nparam
	    dist = ((Memr[parval+ipar-1] - param[ipar])/ scaledist[ipar])** 2 +
		   dist

	# See if this column is the closest to the interpolating value
	# of all in its category

	do icat = 1, ncat {
	    jcat = category[icat]
	    if (collist[jcat] == NULL || dist < mindist[jcat]) {
		collist[jcat] = cp
		mindist[jcat] = dist
	    }
	}

	call sfree (sp)
end




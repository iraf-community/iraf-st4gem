include	<tbset.h>

# COL_SCALE -- Compute scaling factor for parameterized data

procedure col_scale (ncol, colptr, nparam, scaledist)

int	ncol		# i: number of table columns
pointer	colptr[ARB]	# i: array of parameterized columns
int	nparam		# i: number of parameters
real	scaledist[ARB]	# o: scaling factors
#--
double	mean, var
int	icol, ipar, nval
pointer	sp, sum, sumsq, parval, parname

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (sum, nparam, TY_DOUBLE)
	call salloc (sumsq, nparam, TY_DOUBLE)
	call salloc (parval, nparam, TY_REAL)
	call salloc (parname, SZ_COLNAME, TY_CHAR)

	# The scaling factor for a parameter is the standard deviation of
	# the parameter's values

	call amovkd (0.0, Memd[sum], nparam)
	call amovkd (0.0, Memd[sumsq], nparam)

	do icol = 1, ncol {
	    call tbcigt (colptr[icol], TBL_COL_NAME, Memc[parname], SZ_COLNAME)
	    call breakparam (Memc[parname], nparam, nval, Memr[parval])

	    do ipar = 0, nparam-1 {
		Memd[sum+ipar] = Memd[sum+ipar] + Memr[parval+ipar]
		Memd[sumsq+ipar] = Memd[sumsq+ipar] + Memr[parval+ipar] ** 2
	    }
	}

	do ipar = 0, nparam-1 {
	    mean = Memd[sum+ipar] / ncol
	    var = (Memd[sumsq+ipar] - mean * Memd[sum+ipar]) / ncol
	    scaledist[ipar+1] = sqrt (var)
	}
	    
	call sfree (sp)
end

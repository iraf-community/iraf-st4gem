include <tbset.h>

# GETPARCOL -- Get interpolated values from parameterized columns

procedure getparcol (reverse, tp, ncol, colptr, colname, 
		     nparam, param, nrow, wave, value)

bool	reverse		# i: reverse columns before interpolating?
pointer	tp		# i: table descriptor
int	ncol		# i: number of table columns
pointer	colptr[ARB]	# i: table column pointers
char	colname[ARB]	# i: table column name
int	nparam		# i: number of parameters
real	param[ARB]	# i: parameter values
int	nrow		# i: number of rows to read
real	wave[ARB]	# i: wavelength array
real	value[ARB]	# o: column values
#--
int	mcol, icol, ncat
pointer	sp, category, mindist, scaledist, collist

errchk col_category, col_mindist, col_interp

begin
	# Allocate temporary arrays

	mcol = 2 ** nparam

	call smark (sp)
	call salloc (category, mcol, TY_INT)
	call salloc (mindist, mcol, TY_REAL)
	call salloc (collist, mcol, TY_POINTER)
	call salloc (scaledist, nparam, TY_REAL)

	call amovkr (0.0, Memr[mindist], mcol)
	call amovki (NULL, Memi[collist], mcol)

	# Find columns that bracket interpolating value

	call col_scale (ncol, colptr, nparam, Memr[scaledist])

	do icol = 1, ncol {
	    call col_category (colptr[icol], nparam, param, 
			       ncat, Memi[category])

	    call col_mindist (colptr[icol], nparam, param, ncat, 
			      Memi[category], Memr[scaledist], 
			      Memr[mindist], Memi[collist])
	}

	# Perform interpolation

	call col_interp (reverse, tp, nparam, param, Memi[collist], 
			 nrow, wave, value)

	call sfree (sp)
end

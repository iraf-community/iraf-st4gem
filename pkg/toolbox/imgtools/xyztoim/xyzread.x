include <tbset.h>

# xyz_read -- read the data from the input table or ascii file
# The X, Y, and Z data are read from three columns of a table.  The
# values are returned in allocated memory.  The calling routine should
# deallocate this memory when no longer needed.
#
# Phil hodge, 15-Dec-1993  Subroutine created.

procedure xyz_read (intable, xname, yname, zname, ndim,
		x, y, z, nrows)

char	intable[ARB]	# i: name of input table (or ascii file)
char	xname[ARB]	# i: name of column for X values
char	yname[ARB]	# i: name of column for Y values
char	zname[ARB]	# i: name of column for Z values
int	ndim		# i: dimension of image (1 or 2)
pointer x		# o: pointer to X data
pointer y		# o: pointer to Y data
pointer z		# o: pointer to Z data
int	nrows		# o: size of arrays pointed to by x, y, z
#--
pointer tp
pointer cpx, cpy, cpz
double	vx, vy, vz	# x, y, z values for one row
int	allrows		# total number of rows in the table
int	row		# loop index for row number
pointer tbtopn()
int	tbpsta()

begin
	tp = tbtopn (intable, READ_ONLY, NULL)

	call tbcfnd (tp, xname, cpx, 1)
	if (ndim > 1)
	    call tbcfnd (tp, yname, cpy, 1)
	call tbcfnd (tp, zname, cpz, 1)

	if (cpx == NULL || cpz == NULL || (ndim > 1 && cpy == NULL)) {
	    call tbtclo (tp)
	    call error (1, "column not found")
	}

	allrows = tbpsta (tp, TBL_NROWS)

	# Allocate memory for the column data.
	call malloc (x, allrows, TY_DOUBLE)
	if (ndim == 1)
	    call malloc (y, 1, TY_DOUBLE)	# so we can use Memd[y]
	else
	    call malloc (y, allrows, TY_DOUBLE)
	call malloc (z, allrows, TY_DOUBLE)

	# Read the data, ignoring INDEF values, and set nrows.
	vy = 1.d0				# IS_INDEFD(vy) is false
	nrows = 0
	do row = 1, allrows {
	    call tbegtd (tp, cpx, row, vx)
	    if (ndim > 1)
		call tbegtd (tp, cpy, row, vy)
	    call tbegtd (tp, cpz, row, vz)

	    if (!IS_INDEFD(vx) && !IS_INDEFD(vy) && !IS_INDEFD(vz)) {
		Memd[x+nrows] = vx		# zero indexed at this point
		if (ndim > 1)
		    Memd[y+nrows] = vy
		Memd[z+nrows] = vz
		nrows = nrows + 1
	    }
	}

	call tbtclo (tp)
end

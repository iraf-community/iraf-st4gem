include <tbset.h>

# TABCOPY -- Copy the data portion of a fits table

procedure tabcopy (ifile, ofile, title)

char	ifile[ARB]	# i: input file name
char	ofile[ARB]	# i: output file name
bool	title		# i: use first row of table as title?
#--
int	i, ncols, nrows, irow, jrow
pointer	itp, otp, icp, ocp
real	value

int	tbpsta(), tbcnum()
pointer	tbtopn(), open_copy()

errchk	tbtopn, tbegtr, tbeptr

begin
	# Open input and output tables

	itp = tbtopn (ifile, READ_ONLY, NULL)
	otp = open_copy (ofile, itp, title)

	# Create arrays of column pointers for row copy function

	ncols = tbpsta (itp, TBL_NCOLS)
	nrows = tbpsta (itp, TBL_NROWS)

	call salloc (icp, ncols, TY_POINTER)
	call salloc (ocp, ncols, TY_POINTER)

	do i = 1, ncols {
	    Memi[icp+i-1] = tbcnum (itp, i)
	    Memi[ocp+i-1] = tbcnum (otp, i)
	}

	# Copy rows from input to output table

	if (title) {
	    jrow = 2
	    nrows = nrows - 1
	} else {
	    jrow = 1
	}

	do irow = 1, nrows {
	    do i = 1, ncols {
		call tbegtr (itp, Memi[icp+i-1], jrow, value)
		call tbeptr (otp, Memi[ocp+i-1], irow, value)
	    }

	    jrow = jrow + 1
	}

	# Close table

	call tbtclo (itp)
	call tbtclo (otp)

end


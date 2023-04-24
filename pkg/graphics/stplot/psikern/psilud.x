include	<math/iminterp.h>
include	<tbset.h>
include	<psiescape.h>

# Define what columns are needed.
define	RED			1
define	GREEN			2
define	BLUE			3
define	INDEX			4
define	NAME			5
define	NCOLS			5

# Memory management.
define	Blue			Mems[blue+$1]
define	Colname			Memc[colname+(SZ_COLNAME+1)*($1-1)]
define	Colptr			Memi[colptr+$1-1]
define	Green			Mems[green+$1]
define	Iblue			Memr[interp+$1]
define	Igreen			Memr[interp+n+$1]
define	Ired			Memr[interp+2*n+$1]
define	Index			Memi[index+$1-1]
define	Name			Memc[name+(lname+1)*($1)]
define	Null			Memb[null+$1-1]
define	Red			Mems[red+$1]
define	Sblue			Memr[sorted+$1]
define	Sgreen			Memr[sorted+n+$1]
define	Sred			Memr[sorted+2*n+$1]
define	Sx			Memc[sx]
define	Ublue			Memr[unsorted+$1-1]
define	Ugreen			Memr[unsorted+nrows+$1-1]
define	Ured			Memr[unsorted+2*nrows+$1-1]
define	Uname			Memc[uname+(lname+1)*($1-1)]

#---------------------------------------------------------------------------
.help psi_lut_read 4Nov94 source
.ih
NAME
psi_lut_read -- Read in an LUT from a table.
.ih
DESCRIPTION
NOTE: Color indicies start at 0!!!
.endhelp
#---------------------------------------------------------------------------
procedure psi_lut_read (table, red, green, blue, name, n)

char	table[ARB]		# I:  File containing the LUT.
pointer	red			# O:  Pointer to array of red components.
pointer	green			# O:  Pointer to array of green components.
pointer	blue			# O:  Pointer to array of blue components.
pointer	name			# O:  Pointer to array of names.
int	n			# O:  Number of components.

# Declarations
int	ahivi()			# Find maximum of a vector.
pointer	colname			# Name of the columns.
pointer	colptr			# Column descriptors.
int	i			# Generic.
pointer	index			# Index array.
pointer	interp			# Interpolated components.
pointer	lname			# Length of the color names.
int	nrows			# Number of rows in the table.
pointer	null			# Null array.
pointer	sorted			# Sorted components.
pointer	sp			# Stack pointer.
pointer	sx			# Generic string.
pointer	t			# Table to read the colormap from.
pointer	tbcigi()		# Get column parameter.
int	tbpsta()		# Get table parameter.
pointer	tbtopn()		# Open a table
int	type			# Type of table
pointer	uname			# Unsorted color names.
pointer	unsorted		# Unsorted component arrays.

errchk	arbpix
errchk	malloc, mfree
errchk	salloc, sfree, smark
errchk	tbcfnd, tbcgti, tbcgtr, tbcgtt, tbcigi, tbpsta, tbtopn

begin
	call smark (sp)
	call salloc (colname, (SZ_COLNAME+1) * NCOLS, TY_CHAR)
	call salloc (colptr, NCOLS, TY_POINTER)
	call salloc (sx, SZ_LINE, TY_CHAR)
	
	# Open the table.
	t = tbtopn (table, READ_ONLY, NULL)
	nrows = tbpsta (t, TBL_NROWS)
	type = tbpsta (t, TBL_WHTYPE)

	# Set the column names.
	if (type == TBL_TYPE_TEXT) {
	    call strcpy ("c1", Colname(INDEX), SZ_COLNAME)
	    call strcpy ("c2", Colname(RED), SZ_COLNAME)
	    call strcpy ("c3", Colname(GREEN), SZ_COLNAME)
	    call strcpy ("c4", Colname(BLUE), SZ_COLNAME)
	    call strcpy ("c5", Colname(NAME), SZ_COLNAME)
	} else {
	    call strcpy ("red", Colname(RED), SZ_COLNAME)
	    call strcpy ("green", Colname(GREEN), SZ_COLNAME)
	    call strcpy ("blue", Colname(BLUE), SZ_COLNAME)
	    call strcpy ("index", Colname(INDEX), SZ_COLNAME)
	    call strcpy ("name", Colname(NAME), SZ_COLNAME)
	}

	# Find the columns.
	call tbcfnd (t, Colname(1), Colptr(1), NCOLS)

	# Check for column existance.
	do i = 1, NCOLS {
	    if (Colptr(i) == NULL) {
		switch (i) {
		case RED, GREEN, BLUE:
		    call sprintf (Sx, SZ_LINE, "column %s must exist for color component %s")
		    call pargstr (Colname(i))
		    switch(i) {
		    case RED:
			call pargstr ("red")
		    case GREEN:
			call pargstr ("green")
		    case BLUE:
			call pargstr ("blue")
		    }	
		    call error (1, Sx)

		case INDEX:
		    if (type == TBL_TYPE_TEXT) 
			call error (1, "column c1 must exist for the color index in text tables")
		}	
	    }
	}

	# Allocate memory to read from tables.
	call malloc (index, nrows, TY_INT)
	call malloc (null, nrows, TY_BOOL)
	call malloc (unsorted, 3*nrows, TY_REAL)

	# Get the index array.  If non-existant, create one.
	if (Colptr(INDEX) != NULL) {
	    call tbcgti (t, Colptr(INDEX), Index(1), Null(1), 1, nrows)

	    # The index column cannot have any nulls.
	    do i = 1, nrows
		if (Null(i)) {
		    call sprintf (Sx, SZ_LINE,
				  "null value in index column %s at row %d")
		    call pargstr (Colname(INDEX))
		    call pargi (i)
		    call error (1, Sx)
		}
	} else
	    do i = 0, nrows-1
		Index(i) = i

	# Read in RGB.
	call tbcgtr (t, Colptr(RED), Ured(1), Null(1), 1, nrows)
	call tbcgtr (t, Colptr(GREEN), Ugreen(1), Null(1), 1, nrows)
	call tbcgtr (t, Colptr(BLUE), Ublue(1), Null(1), 1, nrows)

	# Read in names.
	if (Colptr(NAME) != NULL) {
	    lname = abs (tbcigi (Colptr(NAME), TBL_COL_DATATYPE))
	    call malloc (uname, (lname+1)*nrows, TY_CHAR)
	    call tbcgtt (t, Colptr(NAME), Uname(1), Null(1), lname, 1, nrows)
	} else {
	    lname = 0
	    call malloc (uname, nrows, TY_CHAR)
	    call amovkc (EOS, Uname(1), nrows)
	}

	# Find maximum of the indicies to decide how many colors are
	# being defined.
	n = ahivi (Index(1), nrows) + 1

	# Allocate arrays for sorted/interpolated LUTs
	call malloc (sorted, 3*n, TY_REAL)
	call amovkr (INDEFR, Memr[sorted], 3*n)
	call malloc (interp, 3*n, TY_REAL)

	# Allocate the resultant arrays.
	call malloc (red, n, TY_SHORT)
	call malloc (green, n, TY_SHORT)
	call malloc (blue, n, TY_SHORT)
	call malloc (name, (lname+1)*n, TY_CHAR)

	# Sort.
	do i = 1, nrows {
	    Sred(Index(i)) = Ured(i)
	    Sgreen(Index(i)) = Ugreen(i)
	    Sblue(Index(i)) = Ublue(i)
	    call strcpy (Uname(i), Name(Index(i)), lname)
	}

	# Interpolate for undefined colors.
	call arbpix (Sred(0), Ired(0), n, II_LINEAR, II_BOUNDARYEXT)
	call arbpix (Sgreen(0), Igreen(0), n, II_LINEAR, II_BOUNDARYEXT)
	call arbpix (Sblue(0), Iblue(0), n, II_LINEAR, II_BOUNDARYEXT)

	# Pack the arrays.
	do i = 0, n-1 {
	    Red(i) = PS_PACKLUT(min (max (Ired(i), 0.), 1.))
	    Green(i) = PS_PACKLUT(min (max (Igreen(i), 0.), 1.))
	    Blue(i) = PS_PACKLUT(min (max (Iblue(i), 0.), 1.))
	}

	# That's all folks.
	call mfree (null, TY_BOOL)
	call mfree (uname, TY_CHAR)
	call mfree (unsorted, TY_REAL)
	call mfree (sorted, TY_REAL)
	call mfree (interp, TY_REAL)
	call mfree (null, TY_BOOL)
	call mfree (index, TY_INT)
	call tbtclo (t)
	call sfree (sp)
end
#---------------------------------------------------------------------------
# End of psi_lut_read
#---------------------------------------------------------------------------

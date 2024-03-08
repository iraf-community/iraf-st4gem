include <imhdr.h>
include <math/gsurfit.h>
include <tbset.h>

# This routine reads an ascii file or a table containing (x, y, z) values,
# fits a surface to z as a function of (x,y), and evaluates the fit at
# (x, y) locations read from a second input table, writing the results of
# the fit to an output table.
#
# Phil Hodge, 20-Dec-1993  Task created.
# RAShaw, 12-Jun-1996	Modified to use table name template package, rather 
#			than fio name template routines.

procedure xyztable()

pointer intable1		# these contain the data to be fit
pointer intable2		# the points at which to evaluate the fit
pointer outtable		# output tables
char	xname[SZ_COLNAME]	# column name for X
char	yname[SZ_COLNAME]	# column name for Y
char	zname[SZ_COLNAME]	# column name for Z
double	gx1, gx2, gy1, gy2	# X and Y at corners (value from cl)
int	xorder, yorder		# number of coefficients to fit
bool	cross_terms		# include cross terms?
pointer sfunction		# function to fit to data
int	func			# function to fit
bool	verbose			# print file names?
bool	coefficients		# print coefficients?
#--
pointer sp
pointer intemplate1	# scratch for list of input tables for data
pointer intemplate2	# scratch for list of input tables for (x,y)
pointer outtemplate	# scratch for list of output tables
pointer sf		# pointer to surface fit struct
pointer x, y, z		# pointers to memory for input data
double	x1, x2, y1, y2	# X and Y at corners of image
double	rms		# RMS deviation of fit from surface
int	ndim		# dimension of fit (1 or 2)
int	nrows		# size of x, y, z arrays
int	i		# loop index for printing coefficients
pointer fit_coeff	# scratch for values returned by dgssave
int	n_coeff		# number of values returned by dgssave
int	dgsgeti()	# for getting n_coeff
double	clgetd()
int	clgeti(), clgwrd()
bool	clgetb()
bool	isblank()
int	inlist1, inlist2, outlist
pointer	tbnopen()		# open a table name template
int	tbnget()		# get next table name from expanded list
int	tbnlen()		# return no. table names in expanded list
int	junk

begin
	call smark (sp)

	call salloc (intemplate1, SZ_FNAME, TY_CHAR)
	call salloc (intemplate2, SZ_FNAME, TY_CHAR)
	call salloc (outtemplate, SZ_FNAME, TY_CHAR)
	call salloc (intable1, SZ_FNAME, TY_CHAR)
	call salloc (intable2, SZ_FNAME, TY_CHAR)
	call salloc (outtable, SZ_FNAME, TY_CHAR)
	call salloc (sfunction, SZ_FNAME, TY_CHAR)

	# Get cl parameters.
	call clgstr ("intable1", Memc[intemplate1], SZ_FNAME)
	if (isblank (Memc[intemplate1]))
	    call error (1, "intable1 is blank")
	call clgstr ("intable2", Memc[intemplate2], SZ_FNAME)
	if (isblank (Memc[intemplate2]))
	    call error (1, "intable2 is blank")
	call clgstr ("outtable", Memc[outtemplate], SZ_FNAME)
	if (isblank (Memc[outtemplate]))
	    call error (1, "outtable is blank")
	call clgstr ("xname", xname, SZ_COLNAME)
	if (isblank (xname))
	    call error (1, "xname must be specified")
	call clgstr ("yname", yname, SZ_COLNAME)
	call clgstr ("zname", zname, SZ_COLNAME)
	if (isblank (zname))
	    call error (1, "zname must be specified")
	xorder = clgeti ("xorder")
	yorder = clgeti ("yorder")
	gx1 = clgetd ("x1")
	gx2 = clgetd ("x2")
	if (isblank (yname) || yorder < 1) {
	    ndim = 1				# 1-D fit
	    yorder = 1
	    gy1 = 0.d0
	    gy2 = 2.d0
	    cross_terms = false
	} else {
	    ndim = 2				# 2-D fit
	    gy1 = clgetd ("y1")
	    gy2 = clgetd ("y2")
	    cross_terms = clgetb ("cross_terms")
	}
	func = clgwrd ("function", Memc[sfunction], SZ_FNAME,
		"|chebyshev|legendre|polynomial")
	verbose = clgetb ("verbose")
	coefficients = clgetb ("coefficients")

	# Convert the first letter to upper case, except for "polynomial".
	if (Memc[sfunction] == 'c')		# Chebyshev
	    Memc[sfunction] = 'C'
	else if (Memc[sfunction] == 'l')	# Legendre
	    Memc[sfunction] = 'L'

	# Open the input and output templates.
	inlist1 = tbnopen (Memc[intemplate1])
	inlist2 = tbnopen (Memc[intemplate2])
	outlist = tbnopen (Memc[outtemplate])

	if (tbnlen (inlist1) != tbnlen (inlist2))
	    call error (1, "Lists of table names must be the same length.")
	if (tbnlen (inlist1) != tbnlen (outlist))
	    call error (1, "Input and output lists must be the same length.")

	while (tbnget (inlist1, Memc[intable1], SZ_FNAME) != EOF) {
	    junk = tbnget (inlist2, Memc[intable2], SZ_FNAME)
	    junk = tbnget (outlist, Memc[outtable], SZ_FNAME)

	    if (verbose) {
		call printf ("%s, %s --> %s; ")		# newline printed later
		    call pargstr (Memc[intable1])
		    call pargstr (Memc[intable2])
		    call pargstr (Memc[outtable])
		call flush (STDOUT)
	    }

	    # Assign values from the cl parameters for this loop, since
	    # these values may be modified for the current table by xyz_limits.
	    x1 = gx1
	    x2 = gx2
	    y1 = gy1
	    y2 = gy2

	    # Read the data from the input table or ascii file.
	    call xyz_read (Memc[intable1], xname, yname, zname, ndim,
		    x, y, z, nrows)

	    # Get values for limits if undefined by user.
	    call xyz_limits (Memd[x], Memd[y], nrows, ndim, x1, x2, y1, y2)

	    # Fit the function to the data.  Also find the RMS deviation.
	    call xyz_fit (sf, Memd[x], Memd[y], Memd[z], nrows, ndim,
		    func, xorder, yorder, cross_terms, x1, x2, y1, y2, rms)
	    if (verbose) {
		call printf (" rms = %.6g\n")
		    call pargd (rms)
		call flush (STDOUT)
	    }

	    if (coefficients) {
		# Write the coefficients.
		n_coeff = dgsgeti (sf, GSNSAVE)
		call salloc (fit_coeff, n_coeff, TY_DOUBLE)
		# Get the values of the coefficients.
		call dgssave (sf, Memd[fit_coeff])
		do i = 0, n_coeff-1 {		# note:  zero indexed
		    call printf ("%.15g\n")
			call pargd (Memd[fit_coeff+i])
		}
		call flush (STDOUT)
	    }

	    # Free memory allocated by xyz_read.
	    call mfree (x, TY_DOUBLE)
	    call mfree (y, TY_DOUBLE)		# OK even if 1-D
	    call mfree (z, TY_DOUBLE)

	    # Read from intable2, write to outtable.
	    call xyz_teval (sf, Memc[intable2], Memc[outtable],
		xname, yname, zname,
		Memc[sfunction], ndim, xorder, yorder, cross_terms, rms)

	    # Free the gsurfit descriptor.
	    call dgsfree (sf)
	}

	call tbnclose (outlist)
	call tbnclose (inlist2)
	call tbnclose (inlist1)
	call sfree (sp)
end

# xyz_teval -- evaluate the fit at each row of the table
# Create the output table as a copy of the input table (intable2).
# At each row of the table, the fit is evaluated, and the value is
# written to the output table.  Parameters are added to the header
# that describe the fit.

procedure xyz_teval (sf, intable, outtable, xname, yname, zname,
		sfunction, ndim, xorder, yorder, cross_terms, rms)

pointer sf		# i: surface fit descriptor
char	intable[ARB]	# i: name of input table containing (x,y)
char	outtable[ARB]	# i: name of output table for (x,y,z)
char	xname[ARB]	# i: column name for X
char	yname[ARB]	# i: column name for Y
char	zname[ARB]	# i: column name for Z
char	sfunction[ARB]	# i: function that was fit to data
int	ndim		# i: dimension of fit
int	xorder, yorder	# i: number of coefficients
bool	cross_terms	# i: were cross terms included in fit?
double	rms		# i: RMS deviation of fit from surface
#--
pointer sp
pointer history			# scratch for history record
pointer itp, otp		# pointers to table descriptors
pointer icpx, icpy		# pointers to column descriptors for intable
pointer ocpz			# pointer to Z column descriptor for outtable
int	row			# row number in input and output tables
int	maxpar			# max number of header parameters
double	x, y			# X and Y at one pixel
double	z			# value of fitted function at (x,y)
pointer tbtopn()
double	dgseval()
int	tbpsta()

begin
	call smark (sp)
	call salloc (history, SZ_FNAME, TY_CHAR)

	itp = tbtopn (intable, READ_ONLY, NULL)

	call tbcfnd (itp, xname, icpx, 1)
	if (icpx == NULL) {
	    call tbtclo (itp)
	    call error (1, "X column not found in intable2")
	}
	if (ndim > 1) {
	    call tbcfnd (itp, yname, icpy, 1)
	    if (icpy == NULL) {
		call tbtclo (itp)
		call error (1, "Y column not found in intable2")
	    }
	}

	# Open output table as a copy of the input.
	otp = tbtopn (outtable, NEW_COPY, itp)

	# Find the Z column if it was in the input table; otherwise, create it.
	call tbcfnd (otp, zname, ocpz, 1)
	if (ocpz == NULL)
	    call tbcdef (otp, ocpz, zname, "", "%20.12g", TY_DOUBLE, 1, 1)

	# Allow room for more header parameters:  xorder, yorder, rmserr,
	# cross_terms, table name, xname, yname, zname, and function.
	if (ndim > 1)
	    maxpar = tbpsta (itp, TBL_MAXPAR) + 9
	else
	    # skip yorder, cross_terms and yname for 1-D
	    maxpar = tbpsta (itp, TBL_MAXPAR) + 6
	call tbpset (otp, TBL_MAXPAR, maxpar)
	call tbtcre (otp)

	# Copy existing header parameters, and write some new ones.
	call tbhcal (itp, otp)			# copy existing keywords
	call tbhadi (otp, "xorder", xorder)
	if (ndim > 1) {
	    call tbhadi (otp, "yorder", yorder)
	    call tbhadb (otp, "xterms", cross_terms)
	}
	call tbhadd (otp, "rmserr", rms)

	call sprintf (Memc[history], SZ_FNAME, "input table name %s")
	    call pargstr (intable)
	call tbhadt (otp, "HISTORY", Memc[history])

	if (ndim > 1) {
	    call sprintf (Memc[history], SZ_FNAME, "X & Y column names %s %s")
		call pargstr (xname)
		call pargstr (yname)
	} else {
	    call sprintf (Memc[history], SZ_FNAME, "X column name %s")
		call pargstr (xname)
	}
	call tbhadt (otp, "HISTORY", Memc[history])
	call sprintf (Memc[history], SZ_FNAME, "Z column name %s")
	    call pargstr (zname)
	call tbhadt (otp, "HISTORY", Memc[history])

	call sprintf (Memc[history], SZ_FNAME, "a %s function was fit")
	    call pargstr (sfunction)
	call tbhadt (otp, "HISTORY", Memc[history])

	# Evaluate the function at each row.  If X or Y is INDEF, set
	# Z to INDEF as well; otherwise, evaluate the function at (X,Y)
	# and save the result in the Z column.
	y = 1.d0				# value for 1-D fit
	do row = 1, tbpsta (itp, TBL_NROWS) {

	    call tbrcpy (itp, otp, row, row)	# copy the row

	    call tbegtd (itp, icpx, row, x)
	    if (ndim > 1)
		call tbegtd (itp, icpy, row, y)

	    if (IS_INDEFD(x) || IS_INDEFD(y)) {
		# Set Z to INDEF as well.
		call tbrudf (otp, ocpz, 1, row)
	    } else {
		# Evaluate the fit at (x,y).
		z = dgseval (sf, x, y)
		call tbeptd (otp, ocpz, row, z)
	    }
	}

	call tbtclo (otp)
	call tbtclo (itp)
	call sfree (sp)
end

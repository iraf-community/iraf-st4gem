include <imhdr.h>
include <math/gsurfit.h>
include <tbset.h>

# This routine reads an ascii file or a table containing (x, y, z) values,
# fits a surface to z as a function of (x,y), and evaluates the fit at
# each point of a uniform rectangular grid, writing the results to an
# output 2-D image.
# If only two columns (x, z) are given, a 1-D image is created.
#
# Phil Hodge, 17-Dec-1993  Task created.
# RAShaw, 12-Jun-1996	Modified to use table name template package, rather 
#			than fio name template routines.

procedure xyztoim()

pointer intable			# list of input table names
pointer output			# list of output images
char	xname[SZ_COLNAME]	# column name for X
char	yname[SZ_COLNAME]	# column name for Y
char	zname[SZ_COLNAME]	# column name for Z
int	nx, ny			# size of image
double	gx1, gx2, gy1, gy2	# X and Y at corners of image (value from cl)
int	xorder, yorder		# number of coefficients to fit
bool	cross_terms		# include cross terms?
pointer sfunction		# function to fit to data
int	func			# function to fit
bool	verbose			# print file names?
bool	coefficients		# print coefficients?
#--
pointer sp
pointer intemplate	# scratch for list of input tables
pointer outtemplate	# scratch for list of output images
pointer im		# pointer to imhdr struct
pointer sf		# pointer to surface fit struct
pointer x, y, z		# pointers to memory for input data
double	x1, x2, y1, y2	# X and Y at corners of image
double	rms		# RMS deviation of fit from surface
int	ndim		# dimension of image (1 or 2)
int	nrows		# size of x, y, z arrays
int	i		# loop index for printing coefficients
pointer fit_coeff	# scratch for values returned by dgssave
int	n_coeff		# number of values returned by dgssave
int	dgsgeti()	# for getting n_coeff
pointer immap()
double	clgetd()
int	clgeti(), clgwrd()
bool	clgetb()
bool	isblank()
int	inlist
pointer	tbnopen()		# open a table name template
int	tbnget()		# get next table name from expanded list
int	tbnlen()		# return no. table names in expanded list
pointer outlist, imtopen()
int	imtlen(), imtgetim(), junk

begin
	call smark (sp)

	call salloc (intemplate, SZ_FNAME, TY_CHAR)
	call salloc (outtemplate, SZ_FNAME, TY_CHAR)
	call salloc (intable, SZ_FNAME, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (sfunction, SZ_FNAME, TY_CHAR)

	# Get cl parameters.
	call clgstr ("intable", Memc[intemplate], SZ_FNAME)
	if (isblank (Memc[intemplate]))
	    call error (1, "no input specified")
	call clgstr ("output", Memc[outtemplate], SZ_FNAME)
	if (isblank (Memc[outtemplate]))
	    call error (1, "no output specified")
	call clgstr ("xname", xname, SZ_COLNAME)
	if (isblank (xname))
	    call error (1, "xname must be specified")
	call clgstr ("yname", yname, SZ_COLNAME)
	call clgstr ("zname", zname, SZ_COLNAME)
	if (isblank (zname))
	    call error (1, "zname must be specified")
	xorder = clgeti ("xorder")
	yorder = clgeti ("yorder")
	nx = clgeti ("nx")
	gx1 = clgetd ("x1")
	gx2 = clgetd ("x2")
	if (isblank (yname) || yorder < 1) {
	    ndim = 1				# create a 1-D image
	    yorder = 1
	    ny = 1
	    gy1 = 0.d0
	    gy2 = 2.d0
	    cross_terms = false
	} else {
	    ndim = 2				# create a 2-D image
	    ny = clgeti ("ny")
	    gy1 = clgetd ("y1")
	    gy2 = clgetd ("y2")
	    cross_terms = clgetb ("cross_terms")
	}
	func = clgwrd ("function", Memc[sfunction], SZ_FNAME,
		"|chebyshev|legendre|polynomial")
	verbose = clgetb ("verbose")
	coefficients = clgetb ("coefficients")

	if (IS_INDEFI(nx) || IS_INDEFI(ny))
	    call error (1, "You must specify nx and ny.")

	# Convert the first letter to upper case, except for "polynomial".
	if (Memc[sfunction] == 'c')		# Chebyshev
	    Memc[sfunction] = 'C'
	else if (Memc[sfunction] == 'l')	# Legendre
	    Memc[sfunction] = 'L'

	# Open the input and output templates.  
	inlist = tbnopen (Memc[intemplate])
	outlist = imtopen (Memc[outtemplate])

	if (tbnlen (inlist) != imtlen (outlist))
	    call error (1, "Input and output lists must be the same length.")

	while (tbnget (inlist, Memc[intable], SZ_FNAME) != EOF) {
	    junk = imtgetim (outlist, Memc[output], SZ_FNAME)

	    if (verbose) {
		call printf ("%s --> %s; ")		# newline printed later
		    call pargstr (Memc[intable])
		    call pargstr (Memc[output])
		call flush (STDOUT)
	    }

	    # Assign values from the cl parameters for this loop, since
	    # these values may be modified for the current table/image.
	    x1 = gx1
	    x2 = gx2
	    y1 = gy1
	    y2 = gy2

	    # Read the data from the input table or ascii file.
	    call xyz_read (Memc[intable], xname, yname, zname, ndim,
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

	    # Create the output image.
	    im = immap (Memc[output], NEW_IMAGE, NULL)
	    IM_NDIM(im) = ndim
	    IM_LEN(im,1) = nx
	    if (ndim > 1)
		IM_LEN(im,2) = ny
	    IM_PIXTYPE(im) = TY_REAL

	    # At each pixel of the image, evaluate the fit and write the
	    # value to the image.
	    call xyz_write (im, sf, x1, x2, y1, y2, nx, ny)

	    # Write info to header.
	    call xyz_header (im, Memc[intable], Memc[sfunction], ndim,
			xorder, yorder, cross_terms, rms, x1, x2, y1, y2,
			xname, yname, zname)

	    # Close the image.
	    call imunmap (im)

	    # Free the gsurfit descriptor.
	    call dgsfree (sf)
	}

	call imtclose (outlist)
	call tbnclose (inlist)
	call sfree (sp)
end

# xyz_write -- evaluate the fit at each pixel of the image
# At each pixel of the image, the fit is evaluated, and the value is
# written to the image at that pixel.  The image is 1-D if ny = 1.
# Parameters are added to the header that describe the mapping from pixel
# coordinates to X and Y.  Datamin & datamax are also set by this routine.

procedure xyz_write (im, sf, x1, x2, y1, y2, nx, ny)

pointer im		# i: pointer to imhdr struct for output image
pointer sf		# i: pointer to gsurfit struct
double	x1, x2		# i: limits on X values
double	y1, y2		# i: limits on Y values
int	nx, ny		# i: size of image (pixels)
#--
pointer ov		# pointer to image data for one line
double	x, y		# X and Y at one pixel
double	dx, dy		# pixel size in X and Y
double	v		# value of fitted function at (x,y)
double	datamin, datamax
int	i, j
pointer impl2r()
double	dgseval()

begin
	# Initial values for min & max; updated in loop.
	datamin = dgseval (sf, x1, y1)
	datamax = datamin

	if (nx > 1)
	    dx = (x2 - x1) / (double(nx) - 1.d0)
	else
	    dx = 1.d0
	if (ny > 1)
	    dy = (y2 - y1) / (double(ny) - 1.d0)
	else
	    dy = 1.d0

	# Evaluate the function, and put the values in the image line by line.
	# Note that the sample (i) and line (j) are zero indexed.
	do j = 0, ny-1 {

	    ov = impl2r (im, j+1)
	    y = y1 + j * dy

	    do i = 0, nx-1 {
		x = x1 + i * dx
		v = dgseval (sf, x, y)		# evaluate the fit at (x,y)
		datamin = min (datamin, v)
		datamax = max (datamax, v)
		Memr[ov+i] = v			# convert to real
	    }
	}

	# Update datamin & datamax header parameters.
        call imputd (im, "i_minpixval", datamin)
        call imputd (im, "i_maxpixval", datamax)
        IM_LIMTIME(im) = IM_MTIME(im) + 1
end

# xyz_header -- add descriptive keywords to image header
# Parameters are added to the header that describe the fit and also
# the mapping from pixel coordinates to X and Y.

procedure xyz_header (im, intable, sfunction, ndim,
		xorder, yorder, cross_terms, rms, x1, x2, y1, y2,
		xname, yname, zname)

pointer im		# i: pointer to imhdr struct for output image
char	intable[ARB]	# i: name of input file
char	sfunction[ARB]	# i: function that was fit to data
int	ndim		# i: dimension of image
int	xorder, yorder	# i: number of coefficients
bool	cross_terms	# i: were cross terms included in fit?
double	rms		# i: RMS deviation of fit from surface
double	x1, x2		# i: limits on X values
double	y1, y2		# i: limits on Y values
char	xname[ARB]	# i: name of column for X values
char	yname[ARB]	# i: name of column for Y values
char	zname[ARB]	# i: name of column for Z values
#--
pointer sp
pointer history		# scratch for history record
double	dx, dy		# pixel size in X and Y
int	nx, ny		# image size

begin
	call smark (sp)
	call salloc (history, SZ_LINE, TY_CHAR)

	nx = IM_LEN(im,1)
	ny = IM_LEN(im,2)
	if (nx > 1)
	    dx = (x2 - x1) / (double(nx) - 1.d0)
	else
	    dx = 1.d0
	if (ny > 1)
	    dy = (y2 - y1) / (double(ny) - 1.d0)
	else
	    dy = 1.d0

	# Add header keywords to give the X, Y coordinates.
	call imaddd (im, "CRPIX1", 1.d0)
	if (ndim > 1)
	    call imaddd (im, "CRPIX2", 1.d0)
	call imaddd (im, "CRVAL1", x1)
	if (ndim > 1)
	    call imaddd (im, "CRVAL2", y1)
	call imaddd (im, "CD1_1", dx)
	if (ndim > 1) {
	    call imaddd (im, "CD1_2", 0.d0)
	    call imaddd (im, "CD2_1", 0.d0)
	    call imaddd (im, "CD2_2", dy)
	}
	call imastr (im, "CTYPE1", "PIXEL")
	if (ndim > 1)
	    call imastr (im, "CTYPE2", "PIXEL")

	call imaddi (im, "XORDER", xorder)
	if (ndim > 1) {
	    call imaddi (im, "YORDER", yorder)
	    call imaddb (im, "XTERMS", cross_terms)
	}
	call imaddd (im, "RMSERR", rms)

	call sprintf (Memc[history], SZ_LINE, "input table name %s")
	    call pargstr (intable)
	call imputh (im, "HISTORY", Memc[history])
	if (ndim > 1) {
	    call sprintf (Memc[history], SZ_LINE, "column names %s %s %s")
		call pargstr (xname)
		call pargstr (yname)
		call pargstr (zname)
	} else {
	    call sprintf (Memc[history], SZ_LINE, "column names %s %s")
		call pargstr (xname)
		call pargstr (zname)
	}
	call imputh (im, "HISTORY", Memc[history])
	call sprintf (Memc[history], SZ_LINE, "a %s function was fit")
	    call pargstr (sfunction)
	call imputh (im, "HISTORY", Memc[history])

	call sfree (sp)
end

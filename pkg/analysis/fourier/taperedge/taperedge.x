include <imhdr.h>
include <error.h>
include <math.h>		# for PI

define	SUBTRACT_NONE	0	# don't subtract anything
define	SUBTRACT_MEAN	1	# subtract mean of entire image
define	SUBTRACT_EDGE	2	# subtract average of edge values
define	SUBTRACT_VALUE	3	# subtract a user-supplied value

define	COS_BELL_TAPER	1	# cosine bell function
define	LINEAR_TAPER	2	# linear taper function

# taperedge -- apply cosine bell function
# This routine applies a cosine bell taper to a list of 1-D or 2-D images.
# The image values will be tapered down from their original values to zero
# at the edges.
#
# The taper affects 'width' pixels at each edge, and the taper is such
# that the result doesn't actually reach zero at the first and last
# pixels.  It would reach zero one pixel beyond the edges of the
# image.  For example, for a 1-D input image filled with a value of one
# using parameters width=5 and subtract="none", the first
# six output values would be:  0.067, 0.25, 0.5, 0.75, 0.933, 1.
#
# This task includes an option to subtract a value before applying the
# cosine bell taper.  The value to be subtracted may be specified as
# "none", "mean", "edge", or a specific numerical value.  Nothing will
# be subtracted if subtract = "none" or "0".  For subtract = "mean",
# the value is the average of the entire input image.  For subtract =
# "edge", the value is the average of the band of 'width' pixels around
# the edge, or for 1-D the average of 'width' pixels at each end.  The
# final alternative is for the user to specify the value to be subtracted.
# After subtracting, the result is the product of the pixel value and
# the taper function.
#
# Phil Hodge, 27-Apr-1993  Task created.

procedure taperedge()

pointer infiles, outfiles	# scratch for names of images (templates)
pointer width			# width of taper zone
pointer sub_str			# for string indicating what to subtract
pointer tfunc_s			# for string giving taper function name
int	subtract		# subtract what?
real	mean			# value to subtract
bool	verbose			# print image names?
#--
pointer sp
pointer input, output		# scratch for specific image names
pointer temp			# scratch for name of temporary image
pointer ilist, olist		# pointers to imt struct
real	widthr			# width of taper zone, possibly in percent
bool	percent			# true if widthr is percent
int	tfunc			# taper function
int	num_output		# number of output image names
int	junk
int	ip
bool	inplace			# true if an image is to be modified in-place
pointer imtopen()
int	imtgetim(), imtlen()
int	ctor()
int	clgwrd()
bool	clgetb()
bool	streq()

begin
	call smark (sp)
	call salloc (infiles, SZ_FNAME, TY_CHAR)
	call salloc (outfiles, SZ_FNAME, TY_CHAR)
	call salloc (width, SZ_FNAME, TY_CHAR)
	call salloc (input, SZ_FNAME, TY_CHAR)
	call salloc (output, SZ_FNAME, TY_CHAR)
	call salloc (sub_str, SZ_FNAME, TY_CHAR)
	call salloc (tfunc_s, SZ_FNAME, TY_CHAR)

	call salloc (temp, SZ_FNAME, TY_CHAR)

	# Get parameters.
	call clgstr ("input", Memc[infiles], SZ_FNAME)
	call clgstr ("output", Memc[outfiles], SZ_FNAME)
	call clgstr ("width", Memc[width], SZ_FNAME)
	call clgstr ("subtract", Memc[sub_str], SZ_FNAME)
	tfunc = clgwrd ("function", Memc[tfunc_s], SZ_FNAME,
		"|cosbell|linear")
	verbose = clgetb ("verbose")

	# Open image templates.
	ilist = imtopen (Memc[infiles])
	olist = imtopen (Memc[outfiles])

	num_output = imtlen (olist)
	if (imtlen (ilist) != num_output && num_output > 0)
	    call error (1, "input and output list are not the same length")

	# Interpret the value to be subtracted.
	mean = 0.
	if (Memc[sub_str] == 'n' || Memc[sub_str] == 'N') {
	    subtract = SUBTRACT_NONE
	} else if (Memc[sub_str] == 'm' || Memc[sub_str] == 'M') {
	    subtract = SUBTRACT_MEAN
	} else if (Memc[sub_str] == 'e' || Memc[sub_str] == 'E') {
	    subtract = SUBTRACT_EDGE
	} else {
	    # Read user-supplied value from string.
	    subtract = SUBTRACT_VALUE
	    ip = 1
	    if (ctor (Memc[sub_str], ip, mean) < 1) {
		call eprintf ("subtract = `%s' is invalid\n")
		    call pargstr (Memc[sub_str])
		call error (1, "")
	    }
	}

	# Read the width from the input string, and check whether the
	# value is a percent or a specific number of pixels.
	call tpr_width (Memc[width], widthr, percent)

	# Do for each image in the input list.
	while (imtgetim (ilist, Memc[input], SZ_FNAME) != EOF) {
	    if (num_output == 0)
		call strcpy (Memc[input], Memc[output], SZ_FNAME)
	    else
		junk = imtgetim (olist, Memc[output], SZ_FNAME)

	    if (streq (Memc[input], Memc[output])) {
		inplace = true
		call mktemp ("taperedge", Memc[temp], SZ_FNAME)
	    } else {
		inplace = false
		call strcpy (Memc[output], Memc[temp], SZ_FNAME)
	    }

	    if (verbose) {
		if (inplace) {
		    call printf ("%s (in-place)\n")
			call pargstr (Memc[input])
		} else {
		    call printf ("%s --> %s\n")
			call pargstr (Memc[input])
			call pargstr (Memc[output])
		}
		call flush (STDOUT)
	    }

	    iferr {
		call bell_taper (Memc[input], Memc[temp], tfunc,
			widthr, percent, subtract, mean, verbose)
	    } then {
		call erract (EA_WARN)
		next
	    }

	    if (inplace) {
		# Delete the input image, and rename the temporary
		# image to the original input name.
		call imdelete (Memc[output])
		call imrename (Memc[temp], Memc[output])
	    }
	}

	call imtclose (olist)
	call imtclose (ilist)

	call sfree (sp)
end

# tpr_width -- read the width from the input string
# The numerical value of width will be read from the width string, and the
# string is searched for "%" or "percent".  

procedure tpr_width (width, widthr, percent)

char	width[ARB]		# i: width of taper zone
real	widthr			# o: width in pixels or percent
bool	percent			# o: true if widthr is in percent
#--
int	ip
int	ctor()
int	strncmp()

begin
	ip = 1
	if (ctor (width, ip, widthr) < 1) {
	    call eprintf ("can't interpret width = `%s'\n")
		call pargstr (width)
	    call error (1, "")
	}

	percent = false				# initial value
	ip = 1
	while (width[ip] != EOS) {
	    if (width[ip] == '%') {
		percent = true
		break
	    } else if (strncmp (width[ip], "percent", 6) == 0) {
		percent = true
		break
	    }
	    ip = ip + 1
	}

	if (percent && widthr > 50.)
	    call error (1, "width greater than 50 percent doesn't make sense")
end

procedure bell_taper (input, output, tfunc,
	widthr, percent, subtract, mean, verbose)

char	input[ARB]	# i: name of input image
char	output[ARB]	# i: name of output image
int	tfunc		# i: taper function
real	widthr		# i: width of taper zone, possibly in percent
bool	percent		# i: true if widthr is percent rather than in pixels
int	subtract	# i: what should we subtract from input image?
real	mean		# io: this is modified if subtract != SUBTRACT_VALUE
bool	verbose		# i: true if we should print info
#--
pointer iim, oim	# imhdr pointers for input & output images
int	xwidth, ywidth	# width of taper in pixels
pointer immap()
errchk	immap, bell_1, bell_2, bell_s1, bell_s2

begin
	iim = immap (input, READ_ONLY, NULL)
	if (IM_NDIM(iim) > 2) {
	    call imunmap (iim)
	    call error (1, "can't handle dimension greater than two")
	}

	if (percent) {
	    xwidth = nint (widthr / 100. * real(IM_LEN(iim,1)))
	    if (2 * xwidth > IM_LEN(iim,1))	# e.g. 50% of odd IM_LEN
		xwidth = xwidth - 1
	    if (IM_NDIM(iim) > 1) {
		ywidth = nint (widthr / 100. * real(IM_LEN(iim,2)))
		if (2 * ywidth > IM_LEN(iim,2))
		    ywidth = ywidth - 1
	    } else {
		ywidth = 1
	    }
	} else {
	    xwidth = nint (widthr)
	    ywidth = nint (widthr)
	}

	if (IM_LEN(iim,1) < 2 * xwidth) {
	    call imunmap (iim)
	    call error (1, "image is smaller than 2 * width")
	}
	if (IM_NDIM(iim) > 1) {
	    if (IM_LEN(iim,2) < 2 * ywidth) {
		call imunmap (iim)
		call error (1, "image is smaller than 2 * width")
	    }
	}

	if (verbose) {
	    if (xwidth == ywidth || IM_NDIM(iim) == 1) {
		call printf ("    width = %d pixels")
		    call pargi (xwidth)
	    } else {
		call printf ("    xwidth = %d, ywidth = %d pixels")
		    call pargi (xwidth)
		    call pargi (ywidth)
	    }
	}
	call flush (STDOUT)

	# Create output image.
	oim = immap (output, NEW_COPY, iim)

	if (subtract == SUBTRACT_NONE) {
	    if (IM_NDIM(iim) == 1)
		call bell_1 (iim, oim, tfunc, xwidth)
	    else if (IM_NDIM(iim) == 2)
		call bell_2 (iim, oim, tfunc, xwidth, ywidth)
	} else {
	    if (IM_NDIM(iim) == 1)
		call bell_s1 (iim, oim, tfunc, xwidth, subtract, mean)
	    else if (IM_NDIM(iim) == 2)
		call bell_s2 (iim, oim, tfunc, xwidth, ywidth, subtract, mean)
	}

	if (verbose && (subtract != SUBTRACT_NONE)) {
	    call printf (";  subtract = %.7g")
		call pargr (mean)
	}
	call printf ("\n")
	call flush (STDOUT)

	call imunmap (oim)
	call imunmap (iim)
end

# For a 1-D image.

procedure bell_1 (iim, oim, tfunc, width)

pointer iim, oim	# i: imhdr pointers for input & output images
int	tfunc		# i: taper function
int	width		# i: width of taper zone
#--
pointer sp
pointer bell		# scratch for cosine bell function
pointer ix, ox		# pointers to image values for input & output
int	npix		# size of image
int	i
pointer imgl1r(), impl1r()

begin
	call smark (sp)
	call salloc (bell, width, TY_REAL)

	# Compute values of taper function.
	call taper_fcn (tfunc, width, Memr[bell])

	npix = IM_LEN(iim,1)

	ix = imgl1r (iim)
	ox = impl1r (oim)

	# First just copy input to output.
	call amovr (Memr[ix], Memr[ox], npix)

	# Now apply the taper.  We'll work forward from the beginning
	# and backward from the end.
	do i = 1, width {
	    Memr[ox+i-1]    = Memr[bell+i-1] * Memr[ix+i-1]
	    Memr[ox+npix-i] = Memr[bell+i-1] * Memr[ix+npix-i]
	}

	call sfree (sp)
end

# For a 2-D image.

procedure bell_2 (iim, oim, tfunc, xwidth, ywidth)

pointer iim, oim	# i: imhdr pointers for input & output images
int	tfunc		# i: taper function
int	xwidth, ywidth	# i: width of taper zone
#--
pointer sp
pointer xbell, ybell	# scratch for cosine bell function
pointer ix, ox		# pointers to image values for input & output
real	w_vert		# cosine bell depending on vertical position in image
real	weight		# value (or product of values in corner) of taper fcn
int	npix1, npix2	# size of image
int	line		# loop index for line of image
int	i
pointer imgl2r(), impl2r()

begin
	call smark (sp)
	call salloc (xbell, xwidth, TY_REAL)
	call salloc (ybell, ywidth, TY_REAL)

	# Compute values of taper function.
	call taper_fcn (tfunc, xwidth, Memr[xbell])
	call taper_fcn (tfunc, ywidth, Memr[ybell])

	npix1 = IM_LEN(iim,1)
	npix2 = IM_LEN(iim,2)

	# Do for each line.
	do line = 1, npix2 {

	    ix = imgl2r (iim, line)
	    ox = impl2r (oim, line)

	    # Get the cosine bell function at this line number.
	    if (line <= ywidth) {
		w_vert = Memr[ybell+line-1]
	    } else if (line > npix2 - ywidth) {
		w_vert = Memr[ybell+npix2-line]
	    } else {
		w_vert = 1.		# we're in the middle region
		call amovr (Memr[ix], Memr[ox], npix1)
	    }

	    # Now apply the taper for the current line.
	    # We'll work forward from the beginning and backward from the end.
	    do i = 1, xwidth {
		weight = Memr[xbell+i-1] * w_vert
		Memr[ox+i-1]     = weight * Memr[ix+i-1]
		Memr[ox+npix1-i] = weight * Memr[ix+npix1-i]
	    }

	    # Near the top and bottom of the image, we need to fill in
	    # the middle portion of the line.
	    if (line <= ywidth || line > (npix2 - ywidth)) {
		weight = w_vert
		do i = xwidth+1, npix1-xwidth
		    Memr[ox+i-1] = weight * Memr[ix+i-1]
	    }
	}

	call sfree (sp)
end

# For a 1-D image.  This version subtracts the mean of the entire image
# or just the edge values.

procedure bell_s1 (iim, oim, tfunc, width, subtract, mean)

pointer iim, oim	# i: imhdr pointers for input & output images
int	tfunc		# i: taper function
int	width		# i: width of taper zone
int	subtract	# i: what should we subtract from input image?
real	mean		# io: this is input if subtract = SUBTRACT_VALUE
#--
pointer sp
pointer bell		# scratch for cosine bell function
pointer ix, ox		# pointers to image values for input & output
real	sum		# for subtracting mean
int	npix		# size of image
int	i
pointer imgl1r(), impl1r()

begin
	call smark (sp)
	call salloc (bell, width, TY_REAL)

	# Compute values of taper function.
	call taper_fcn (tfunc, width, Memr[bell])

	npix = IM_LEN(iim,1)

	ix = imgl1r (iim)
	ox = impl1r (oim)

	# Find the mean.
	sum = 0.
	if (subtract == SUBTRACT_MEAN) {
	    do i = 0, npix-1			# note:  zero indexed
		sum = sum + Memr[ix+i]
	    mean = sum / npix
	} else if (subtract == SUBTRACT_EDGE) {
	    do i = 1, width {
		sum = sum + Memr[ix+i-1]
		sum = sum + Memr[ix+npix-i]
	    }
	    if (width > 0)
		mean = sum / (2. * width)
	    else
		mean = 0.
	} else if (subtract == SUBTRACT_VALUE) {
	    # We already have the value to subtract.
	    ;
	} else {
	    mean = 0.
	}

	# Subtract the mean, writing to output.
	call asubkr (Memr[ix], mean, Memr[ox], npix)

	# Now apply the taper.  We'll work forward from the beginning
	# and backward from the end.
	do i = 1, width {
	    Memr[ox+i-1]    = Memr[bell+i-1] * Memr[ox+i-1]
	    Memr[ox+npix-i] = Memr[bell+i-1] * Memr[ox+npix-i]
	}

	call sfree (sp)
end

# For a 2-D image.  This version subtracts the mean of the entire image
# or just the edge values.

procedure bell_s2 (iim, oim, tfunc, xwidth, ywidth, subtract, mean)

pointer iim, oim	# i: imhdr pointers for input & output images
int	tfunc		# i: taper function
int	xwidth, ywidth	# i: width of taper zone
int	subtract	# i: what should we subtract from input image?
real	mean		# io: this is input if subtract = SUBTRACT_VALUE
#--
pointer sp
pointer xbell, ybell	# scratch for cosine bell function
pointer ix, ox		# pointers to image values for input & output
real	sum		# for subtracting mean
real	num_edge	# for subtracting mean of pixel values around edge
real	w_vert		# cosine bell depending on vertical position in image
real	weight		# value (or product of values in corner) of taper fcn
int	npix1, npix2	# size of image
int	line		# loop index for line of image
int	i
pointer imgl2r(), impl2r()

begin
	call smark (sp)
	call salloc (xbell, xwidth, TY_REAL)
	call salloc (ybell, ywidth, TY_REAL)

	# Compute values of taper function.
	call taper_fcn (tfunc, xwidth, Memr[xbell])
	call taper_fcn (tfunc, ywidth, Memr[ybell])

	npix1 = IM_LEN(iim,1)
	npix2 = IM_LEN(iim,2)

	# Find the mean.
	sum = 0.

	if (subtract == SUBTRACT_MEAN) {

	    do line = 1, npix2 {
		ix = imgl2r (iim, line)

		do i = 0, npix1-1		# note:  zero indexed
		    sum = sum + Memr[ix+i]
	    }
	    mean = sum / (npix1 * npix2)

	} else if (subtract == SUBTRACT_EDGE) {

	    if (xwidth > 0 || ywidth > 0) {
		do line = 1, npix2 {
		    ix = imgl2r (iim, line)

		    if (line <= ywidth || line > npix2 - ywidth) {
			do i = 0, npix1-1 {		# note:  zero indexed
			    sum = sum + Memr[ix+i]
			}
		    } else {
			do i = 1, xwidth {
			    sum = sum + Memr[ix+i-1]
			    sum = sum + Memr[ix+npix1-i]
			}
		    }
		}
		num_edge = 2 * (ywidth * npix1 + xwidth * npix2) -
				4 * xwidth * ywidth
		mean = sum / real (num_edge)
	    } else {
		mean = 0.
	    }
	} else if (subtract == SUBTRACT_VALUE) {
	    # We already have the value to subtract.
	    ;
	} else {
	    mean = 0.
	}

	# Do for each line.
	do line = 1, npix2 {

	    ix = imgl2r (iim, line)
	    ox = impl2r (oim, line)

	    # Subtract the mean, writing to output.
	    call asubkr (Memr[ix], mean, Memr[ox], npix1)

	    # Get the weight for this line number.
	    if (line <= ywidth)
		w_vert = Memr[ybell+line-1]
	    else if (line > npix2 - ywidth)
		w_vert = Memr[ybell+npix2-line]
	    else
		w_vert = 1.		# we're in the middle region

	    # Now apply the taper for the current line.
	    # Work forward from the beginning and backward from the end.
	    do i = 1, xwidth {
		weight = Memr[xbell+i-1] * w_vert
		Memr[ox+i-1]     = weight * Memr[ox+i-1]
		Memr[ox+npix1-i] = weight * Memr[ox+npix1-i]
	    }

	    # Near the top and bottom of the image, we also need to fill in
	    # the middle portion of the line.
	    if (line <= ywidth || line > (npix2 - ywidth)) {
		weight = w_vert
		do i = xwidth+1, npix1-xwidth
		    Memr[ox+i-1] = weight * Memr[ox+i-1]
	    }
	}

	call sfree (sp)
end

# taper_fcn -- compute the taper function
# This routine computes the values of the taper function.  The first
# element of the output array is the smallest value (near zero), and the
# last element is the largest value (near one).

procedure taper_fcn (tfunc, width, bell)

int	tfunc		# i: indicates which function to use
int	width		# i: width of taper zone
real	bell[width]	# o: the values of the taper function
#--
real	arg		# function argument
int	i

begin
	if (tfunc == COS_BELL_TAPER) {
	    do i = 1, width {
		arg = real(i) * PI / real (width + 1)
		bell[i] = (1. - cos (arg)) / 2.
	    }
	} else if (tfunc == LINEAR_TAPER) {
	    do i = 1, width
		bell[i] = real(i) / real (width + 1)
	}
end

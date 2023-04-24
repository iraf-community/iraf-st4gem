include	<imio.h>
include	<imhdr.h>

define	MAX_IMAGE	128

# The procedures in this file handle the image variables used in the 
# expression. First, the task calls init_imvar() which opens the images
# and initializes the line counter. Then the task calls next_imvar() 
# which reads the next line from each image and increments the line 
# counter. Then vex_eval() calls get_imvar() to get a line from an
# image as it needs to. Finally, all the images are closed and memory
# released by free_imvar().
#
# B.Simon	18-May-90	Original
# B.Simon	24-Apr-91	Revised to handle multiple data types
# B.Simon	13-Jul-93	Revised to handle multi-group images
# B.Simon	26-Jan-96	Don't open images not mentioned in expression
# B.Simon	17-Oct-97	Catch errors that occur during error handling
# B.Simon	23-Oct-97	Get rid of err_imvar
# B.Simon	06-Jul-98	Check dimensions of image opened read/write
# B.Simon	14-Sep-98	Changes EOF exit code in next_imvar

# FIND_IMVAR -- Find which images variables are mentioned in the expression

procedure find_imvar (equals, found, nfound)

char	equals[ARB]	# i: Expression to be evaluated
int	found[ARB]	# o: The array of images found in the expression
int	nfound		# i: Size of the array
#--
int	fd, ic, jc, nc, ival
pointer	sp, line

int	open(), getline(), strsearch(), ctoi()

begin
	# Allocate temporary array to hold line

	call smark (sp)
	call salloc (line, SZ_LINE, TY_CHAR)

	# Is line a file?

	for (ic = 1; equals[ic] != EOS && equals[ic] <= ' '; ic = ic + 1)
	    ;

	if (equals[ic] == '@') {
	    fd = open (equals[ic+1], READ_ONLY, TEXT_FILE)
	    nc = getline (fd, Memc[line])
	} else {
	    fd = 0
	    call strcpy (equals, Memc[line], SZ_LINE)
	}

	# Initialize array elements to no

	call amovki (NO, found, nfound)

	repeat {
	    ic = 1

	    repeat {
		# Find the next image variable in the expression

		jc = strsearch (Memc[line+ic-1], "im")
		if (jc == 0)
		    break

		# Set the corresponding array element to yes

		ic = ic + jc - 1
		nc = ctoi (Memc[line], ic, ival)
		if (ival > 0 && ival < nfound)
		    found[ival+1] = YES
	    }

	    # Termination criterion for outer loop

	    if (fd == 0) {
		break
	    } else if (getline (fd, Memc[line]) == EOF) {
		call close (fd)
		break
	    }
	}

	call sfree (sp)
end

# FREE_IMVAR -- Free the image variable arrays

procedure free_imvar ()

#--
include	"imvar.com"

int	index
errchk	gf_unmap

begin
	# Check to see if routine has already been called

	if (imptr == NULL)
	    return

	# Close each input image

	do index = 0, nimage  {
	    if (Memi[imptr+index] != NULL)
		call gf_unmap (Memi[imptr+index])	## gf_unmap
	    Memi[imptr+index] = NULL
	}
	
	call mfree (imptr, TY_INT)
	call mfree (dataptr, TY_INT)
	call mfree (typeptr, TY_INT)
	call mfree (gindex, TY_INT)
	call mfree (gcount, TY_INT)
	call mfree (oldline, TY_LONG)
	call mfree (newline, TY_LONG)
	call mfree (outroot, TY_CHAR)
	call mfree (outsect, TY_CHAR)

	imptr = NULL

end

# GET_IMVAR -- Get a line from the named image

procedure get_imvar (stack, name)

pointer	stack		# u: Expression stack pointer
char	name[ARB]	# i: Image name
#--
include	"imvar.com"

bool	found
int	ic, index, junk, ipix
pointer	sp, errmsg, buffer

string	badname  "Illegal variable name (%s)"

int	ctoi()
pointer	stk_alloc()

begin
	call smark (sp)
	call salloc (errmsg, SZ_LINE, TY_CHAR)


	# Determine the kind of variable from its name

	found = false
	if (name[1] == 'i' && name[2] == 'm') { 
	    # Variable is an image line

	    ic = 3
	    junk = ctoi (name, ic, index)

	    if (index > 0 && index <= nimage) {
		if (Memi[dataptr+index] != NULL) {
		    found = true
		    switch (Memi[typeptr+index]) {
		    case TY_INT, TY_LONG:
			buffer = stk_alloc (stack, npix, TY_INT)
			call amovi (Memi[Memi[dataptr+index]], 
				    Memi[buffer], npix)
		    case TY_REAL:
			buffer = stk_alloc (stack, npix, TY_REAL)
			call amovr (Memr[Memi[dataptr+index]], 
				    Memr[buffer], npix)
		    case TY_DOUBLE:
			buffer = stk_alloc (stack, npix, TY_DOUBLE)
			call amovd (Memd[Memi[dataptr+index]], 
				    Memd[buffer], npix)
		    }
		}
	    } 

	} else if (name[2] == EOS) {
	    # Variable is a line number

	    found = true
	    switch (name[1]) {
	    case 'x':	    
		buffer = stk_alloc (stack, npix, TY_INT)
		do ipix = 1, npix
		    Memi[buffer+ipix-1] = ipix

	    case 'y':
		buffer = stk_alloc (stack, npix, TY_INT)
		do ipix = 1, npix
		    Memi[buffer+ipix-1] = Meml[oldline+1]

	    case 'z':
		buffer = stk_alloc (stack, npix, TY_INT)
		do ipix = 1, npix
		    Memi[buffer+ipix-1] = Meml[oldline+2]

	    default:
		found = false
	    }
	}

	if (! found) {
	    call sprintf (Memc[errmsg], SZ_LINE, badname)
		call pargstr (name)
	    call error (1, Memc[errmsg])
	}

	call sfree (sp)

end

# GRP_IMVAR -- Read next group from image

int procedure grp_imvar (index)

int	index		# i: Index to group number
#--
include	"imvar.com"

real	datamin, datamax
pointer	sp, image

pointer	gf_map()

begin
	# Check group index against total number of groups

	if (Memi[gcount+index] <= 1)
	    return (NO)

	# Allocate memory for temporary string

	call smark (sp)
	call salloc (image, SZ_FNAME, TY_CHAR)

	# Increment group index, decrement group count and set line number to 1

	Memi[gindex+index] = Memi[gindex+index] + 1
	Memi[gcount+index] = Memi[gcount+index] - 1
	call amovkl (long(1), Meml[newline], IM_MAXDIM)

	# If this is an input image, update the group parameter block
	# in the existing descriptor. If it is the output image create
	# the new image name, close the old image, and open the new

	if (index != 0) {
	    call gf_opengr (Memi[imptr+index], Memi[gindex+index], 
			    datamin, datamax, 0)
	} else {
	    call sprintf (Memc[image], SZ_FNAME, "%s[%d]%s")
	    call pargstr (Memc[outroot])
	    call pargi (Memi[gindex])
	    call pargstr (Memc[outsect])

	    call gf_unmap (Memi[imptr])	## gf_unmap
	    Memi[imptr] = gf_map (Memc[image], NEW_COPY, Memi[imptr+1])
	}

	call sfree (sp)
	return (YES)
end

# INIT_IMVAR -- Initialize image variable arrays

int procedure init_imvar (input, output, equals, pixtype)

char	input[ARB]	# i: Input image name template
char	output[ARB]	# i: Output image name
char	equals[ARB]	# i: Expression to be evaluated
int	pixtype		# u: Output image pixel type
#--
include	"imvar.com"

int	overwrite, jimage, ndimen, idim, pass, totline, status
pointer	sp, image, root, section, found, list, buffer

string	badgroup "Images do not have same number of groups"
string	badsize  "Image dimensions are not the same"

bool	streq()
int	imtgetim()
pointer	imtopen(), gf_map(), impnli(), impnlr(), impnld()

begin

	# Allocate memory for strings and arrays

	call smark (sp)
	call salloc (image, SZ_FNAME, TY_CHAR)
	call salloc (root, SZ_FNAME, TY_CHAR)
	call salloc (section, SZ_FNAME, TY_CHAR)
	call salloc (found, MAX_IMAGE, TY_INT)

	# Allocate pointer arrays

	call malloc (imptr, MAX_IMAGE, TY_INT)
	call malloc (dataptr, MAX_IMAGE, TY_INT)
	call malloc (typeptr, MAX_IMAGE, TY_INT)
	call malloc (gindex, MAX_IMAGE, TY_INT)
	call malloc (gcount, MAX_IMAGE, TY_INT)
	call malloc (oldline, IM_MAXDIM, TY_LONG)
	call malloc (newline, IM_MAXDIM, TY_LONG)
	call malloc (outroot, SZ_FNAME, TY_CHAR)
	call malloc (outsect, SZ_FNAME, TY_CHAR)

	# Determine which image variables are used in the expression

	call find_imvar (equals, Memi[found], MAX_IMAGE)

	# Process each input image

	pass = 1
	nline = 0
	overwrite = NO

	list = imtopen (input)
	call imgcluster (output, Memc[outroot], SZ_FNAME)

	repeat {
	    jimage = 0
	    nimage = 0

	    while (imtgetim (list, Memc[image], SZ_FNAME) != EOF) {

		nimage = nimage + 1
		if (nimage >= MAX_IMAGE)
		    call error (1, "Maximum number of images exceeded")

		# Check to see if image is in the expression, if not, 
		# set its pointers to NULL and process next image

		# Fix: Make two passes over list of images and only
		# check to see if they are in expression on the first
		# pass. This handles the case when we have an expression
		# with no images (a constant expression) but need to open
		# an image to determine the size of the output image

		if (pass == 1 && Memi[found+nimage] == NO) {
		    Memi[imptr+nimage] = NULL
		    Memi[dataptr+nimage] = NULL
		    Memi[typeptr+nimage] = NULL
		    next
		}

		# Open input image

		Memi[imptr+nimage] = gf_map (Memc[image], READ_ONLY, 0)
		Memi[dataptr+nimage] = NULL

		# Get group index and number of groups

		call tp_parse (Memc[image], 0, Memc[root], Memc[section],
			       Memi[gindex+nimage], Memi[gcount+nimage])

		if (streq (Memc[root], Memc[outroot]))
		    overwrite = YES

		# Set image type

		switch (IM_PIXTYPE(Memi[imptr+nimage])) {
		case TY_SHORT,TY_INT,TY_LONG,TY_USHORT:
		    Memi[typeptr+nimage] = TY_INT
		case TY_REAL:
		    Memi[typeptr+nimage] = TY_REAL
		case TY_DOUBLE:
		    Memi[typeptr+nimage] = TY_DOUBLE
		default:
		    Memi[typeptr+nimage] = TY_REAL
		}

		# Check that dimensions and number of groups are equal

		if (jimage == 0) {
		    jimage = nimage

		    ndimen = IM_NDIM(Memi[imptr+jimage])
		    npix = IM_LEN(Memi[imptr+jimage],1)
		    call amovl (IM_LEN(Memi[imptr+jimage],1), 
				Meml[oldline], IM_MAXDIM)

		} else {

		    if (ndimen != IM_NDIM(Memi[imptr+nimage]))
			call error (1, badsize)

		    do idim = 1, ndimen {
			if (IM_LEN(Memi[imptr+nimage],idim) !=
			    Meml[oldline+idim-1])
			    call error (1, badsize)
		    }

		    if (Memi[gcount+nimage] != Memi[gcount+jimage])
			call error (1, badgroup)
		}
	    }

	    if (nimage == 0)
		call error (1, "Input image not found")

	    pass = pass + 1
	    call imtrew (list)

	} until (jimage > 0)

	call imtclose (list)

	# Parse output image name to get group index and number of groups

	call tp_parse (output, Memi[gcount+jimage], Memc[outroot], 
		       Memc[outsect], Memi[gindex], Memi[gcount])

	# Open output image

	if (overwrite == YES) {
	    Memi[imptr] = gf_map (output, READ_WRITE, NULL)

	    if (ndimen != IM_NDIM(Memi[imptr]))
		call error (1, badsize)

	    do idim = 1, ndimen {
		if (IM_LEN(Memi[imptr],idim) !=
		    Meml[oldline+idim-1])
		    call error (1, badsize)
	    }

	    if (Memi[gcount] != Memi[gcount+jimage])
		call error (1, badgroup)

	} else {
	    if (Memi[gcount] == 1) {
		call strcpy (output, Memc[image], SZ_FNAME)

	    } else {
		call sprintf (Memc[image], SZ_FNAME, "%s[%d/%d]%s")
		call pargstr (Memc[outroot])
		call pargi (Memi[gindex])
		call pargi (Memi[gcount])
		call pargstr (Memc[outsect])
	    }

	    Memi[imptr] = gf_map (Memc[image], NEW_COPY, Memi[imptr+jimage])
	}

	Memi[dataptr] = NULL

	if (pixtype > 0 && overwrite == NO) {
	    IM_PIXTYPE(Memi[imptr]) = pixtype
	} else {
	    pixtype = IM_PIXTYPE(Memi[imptr+jimage])
	}

	# Compute total number of lines
	totline = Memi[gcount+jimage]
	do idim = 2, ndimen
	   totline = totline * Meml[oldline+idim-1]

	# Set the pixel type and write a phony line to the image 
	# to force the creation of the data file

	call amovkl (long(1), Meml[oldline], IM_MAXDIM)

	switch (pixtype) {
	case TY_SHORT,TY_INT,TY_LONG,TY_USHORT:
	    Memi[typeptr] = TY_INT
	    status = impnli (Memi[imptr], buffer, Meml[oldline])
	    call amovki (0, Memi[buffer], npix)
	case TY_REAL:
	    Memi[typeptr] = TY_REAL
	    status = impnlr (Memi[imptr], buffer, Meml[oldline])
	    call amovkr (0.0, Memr[buffer], npix)
	case TY_DOUBLE:
	    Memi[typeptr] = TY_DOUBLE
	    status = impnld (Memi[imptr], buffer, Meml[oldline])
	    call amovkd (0.0d0, Memd[buffer], npix)
	default:
	    Memi[typeptr] = TY_REAL
	    status = impnlr (Memi[imptr], buffer, Meml[oldline])
	    call amovkr (0.0, Memr[buffer], npix)
	}

	call amovkl (long(1), Meml[oldline], IM_MAXDIM)
	call amovkl (long(1), Meml[newline], IM_MAXDIM)
	return (totline)
	
end

# NEXT_IMVAR -- Get next line from the images

int procedure next_imvar (buffer, npixel)

pointer	buffer		# o: Pointer to output buffer
int	npixel		# o: Length of output line
#--
include	"imvar.com"

bool	done
int	index, status

int	imgnli(), impnli(), imgnlr(), impnlr(), imgnld(), impnld()
int	grp_imvar()

begin
	# Set length of output buffer

	npixel = npix

	# Set the previous line to the current line

	call amovl (Meml[newline], Meml[oldline], IM_MAXDIM)

	# Read the input images. Because the input procedure
	# increments the current line number, the current line
	# number must be copied from the previous line number
	# at each pass

	done = false
	do index = 1, nimage {
	    call amovl (Meml[oldline], Meml[newline], IM_MAXDIM)
	    repeat {
		switch (Memi[typeptr+index]) {
		case NULL:  
		    ;		# image was not opened
		case TY_INT, TY_LONG:
		    status = imgnli (Memi[imptr+index], Memi[dataptr+index], 
				     Meml[newline])
		case TY_REAL:
		    status = imgnlr (Memi[imptr+index], Memi[dataptr+index], 
				     Meml[newline])
		case TY_DOUBLE:
		    status = imgnld (Memi[imptr+index], Memi[dataptr+index], 
				     Meml[newline])
		}

		if (status == EOF) {
		    if (grp_imvar (index) == NO) {
			done = true
			status = OK
		    }
		}

	    } until (status != EOF)
	}

	# Get the data buffer for the output image. The output
	# procedure also increments the current line number

	call amovl (Meml[oldline], Meml[newline], IM_MAXDIM)
	repeat {
	    switch (Memi[typeptr]) {
	    case TY_INT, TY_LONG:
		status = impnli (Memi[imptr], buffer, Meml[newline])
	    case TY_REAL:
		status = impnlr (Memi[imptr], buffer, Meml[newline])
	    case TY_DOUBLE:
		status = impnld (Memi[imptr], buffer, Meml[newline])
	    }

	    if (status == EOF) {
		if (grp_imvar (0) == NO) {
		    done = true
		    status = OK
		}
	    }

	} until (status != EOF)

	if (done) {
	    return (EOF)
	} else {
	    nline = nline + 1
	    return (nline)
	}
end

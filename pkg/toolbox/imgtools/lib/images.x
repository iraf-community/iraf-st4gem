include <imio.h>
include <imhdr.h>
include	<ctype.h>

# These routines perform i/o to multigroup images and bad pixel lists,
# disguising the differences between the two so that higher level 
# routines need not worry about the differences

define	LEN_IMGSTRUCT	6

define	IMG_DESCRIP	Memi[$1]		# image descriptor
define	IMG_NAMPTR	Memi[$1+1]		# pointer to image name
define	IMG_SECTPTR	Memi[$1+2]		# pointer to image section
define	IMG_PIXLIST	Memi[$1+3]		# is image a bad pixel list?
define	IMG_MODE	Memi[$1+4]		# image access mode
define	IMG_GCOUNT	Memi[$1+5]		# total number of groups

define	IMG_NAME	Memc[IMG_NAMPTR($1)]	# image name
define	IMG_SECTION	Memc[IMG_SECTPTR($1)]	# image section

define	NFIELD		4			# number of fields in name
define	TYP_KWORD	"PIXTYPE"		# keyword name of pixel type

# CLS_IMAGE -- Close image and free associated structure

procedure cls_image (img)

pointer	img		# i: image descriptor
#--
pointer	im

errchk	imunmap

begin
	# Free image descriptor

	im = IMG_DESCRIP(img)

	call mfree (IMG_NAMPTR(img), TY_CHAR)
	call mfree (IMG_SECTPTR(img), TY_CHAR)
	call mfree (img, TY_INT)

	# Close image

	if (im != NULL)
	    call imunmap (im)

end

# CPY_IMAGE -- Create a new image which is a copy of an old one

procedure cpy_image (name, ndim, oldimg, gindex, newimg)

char	name[ARB]	# i: name of new image
int	ndim		# i: dimensionality of new image (0 means same as old)
pointer	oldimg		# i: old descriptor
int	gindex		# o: index to current group
pointer	newimg		# o: new descriptor
#--
int	idim, dim1, dim2, gcount1, gcount2, type
pointer	sp, size, fullname, errmsg, im2

string	nogroup  "Output image does not support group format (%s)"
string	nomatch  "Number of groups in output does not match input (%s)"

bool	chkexten()
pointer	immap()

begin
	call smark (sp)
	call salloc (size, IM_MAXDIM, TY_LONG)
	call salloc (fullname, SZ_FNAME, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	# Allocate image structure

	call malloc (newimg, LEN_IMGSTRUCT, TY_INT)
	call malloc (IMG_NAMPTR(newimg), SZ_FNAME, TY_CHAR)
	call malloc (IMG_SECTPTR(newimg), SZ_FNAME, TY_CHAR)

	# Get information from image fields

	call parse_image (name, gindex, gcount2, 
			 IMG_NAME(newimg), IMG_SECTION(newimg), SZ_FNAME)

	# Set the number of groups and the dimensions of the output image

	call par_image (oldimg, type, dim1, gcount1, Meml[size])

	if (! chkexten (name, "??h") || chkexten (name, "imh")) {
	    if (gcount2 == 0) {
		gcount2 = 1
	    } else if (gcount2 > 1) {
		call sprintf (Memc[errmsg], SZ_LINE, nogroup)
		call pargstr (name)
		call error (1, Memc[errmsg])
	    }
	}

	if (gcount2 == 0) {
	    if (ndim == 0 || ndim == dim1) {
		dim2 = dim1
		gcount2 = gcount1
	    } else if (ndim < dim1) {
		dim2 = ndim
		gcount2 = gcount1
		do idim = dim2+1, dim1 {
		    gcount2 = gcount2 * Meml[size+idim-1]
		    Meml[size+idim-1] = 1
		}
	    } else {
		dim2 = min (ndim, IM_MAXDIM)
		gcount2 = 1
		Meml[size+dim1] = gcount1
	    }

	} else if (gcount1 == gcount2) {
	    dim2 = dim1

	} else if (gcount2 == 1) {
	    if (ndim == 0) {
		dim2 = dim1 + 1
		Meml[size+dim1] = gcount1
	    } else if (ndim > dim1) {
		dim2 = min (ndim, IM_MAXDIM)
		Meml[size+dim1] = gcount1
	    } else {
		call sprintf (Memc[errmsg], SZ_LINE, nomatch)
		call pargstr (name)
		call error (1, Memc[errmsg])
	    }

	} else {
	    call sprintf (Memc[errmsg], SZ_LINE, nomatch)
	    call pargstr (name)
	    call error (1, Memc[errmsg])
	}

	# Rebuild image name from parts

	if (gcount2 == 1) {
	    call strcpy (IMG_NAME(newimg), Memc[fullname], SZ_FNAME)
	    call strcat (IMG_SECTION(newimg), Memc[fullname], SZ_FNAME)

	} else {
	    call sprintf (Memc[fullname], SZ_FNAME, "%s[%d/%d]%s")
	    call pargstr (IMG_NAME(newimg))
	    call pargi (gindex)
	    call pargi (gcount2)
	    call pargstr (IMG_SECTION(newimg))
	}

	# Check to see if new image is a pixel list
	# if so, open it as a new file, if not,
	# open it as a new copy of the old image

	if (chkexten (name, "pl")) {
	    IMG_PIXLIST(newimg) = YES
	    IMG_MODE(newimg) = NEW_FILE

	    im2 = immap (Memc[fullname], NEW_FILE, NULL)
	    call pl_ssize (IM_PL(im2), dim2, Meml[size], 0)
	    call imaddi (im2, TYP_KWORD, type)

	} else {
	    IMG_PIXLIST(newimg) = NO
	    IMG_MODE(newimg) = NEW_COPY

	    im2 = immap (Memc[fullname], NEW_COPY, IMG_DESCRIP(oldimg))
	    IM_PIXTYPE(im2) = type
	}

	# Set image dimensions


	IM_NDIM(im2) = dim2
	IM_NPHYSDIM(im2) = dim2

	do idim = 1, IM_MAXDIM {
	    IM_LEN(im2, idim) = Meml[size+idim-1]
	    IM_SVLEN(im2,idim) = Meml[size+idim-1]
	}

	# Fill in remaining fields of image descriptor

	IMG_DESCRIP(newimg) = im2
	IMG_GCOUNT(newimg) = gcount2

	call sfree (sp)
end

# DEL_IMAGE -- Delete image and free associated structure

procedure del_image (img)

pointer	img		# i: image descriptor
#--
pointer	im

errchk	imunmap, imdelete

begin
	# Delete image

	im = IMG_DESCRIP(img)

	if (im != NULL) {
	    call imunmap (im)
	    call imdelete (IMG_NAME(img))
	}

	# Free image descriptor

	call mfree (IMG_NAMPTR(img), TY_CHAR)
	call mfree (IMG_SECTPTR(img), TY_CHAR)
	call mfree (img, TY_INT)

end

# NEW_IMAGE -- Create a new image which is a copy of a table

procedure new_image (name, ndim, tab, gindex, img)

char	name[ARB]	# i: name of new image
int	ndim		# i: dimensionality of new image (0 means same as old)
pointer	tab		# i: table descriptor
int	gindex		# o: index to current group
pointer	img		# o: image descriptor
#--
int	idim, dim1, dim2, gcount1, gcount2, type
pointer	sp, size, fullname, errmsg, im2

string	nogroup  "Output image does not support group format (%s)"
string	nomatch  "Number of groups in output does not match input (%s)"

bool	chkexten()
pointer	immap()

begin
	call smark (sp)
	call salloc (size, IM_MAXDIM, TY_LONG)
	call salloc (fullname, SZ_FNAME, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	# Allocate image structure

	call malloc (img, LEN_IMGSTRUCT, TY_INT)
	call malloc (IMG_NAMPTR(img), SZ_FNAME, TY_CHAR)
	call malloc (IMG_SECTPTR(img), SZ_FNAME, TY_CHAR)

	# Get information from image fields

	call parse_image (name, gindex, gcount2, 
			 IMG_NAME(img), IMG_SECTION(img), SZ_FNAME)

	# Set the number of groups and the dimensions of the output image

	call par_table (tab, type, dim1, gcount1, Meml[size])

	if (! chkexten (name, "??h") || chkexten (name, "imh")) {
	    if (gcount2 == 0) {
		gcount2 = 1
	    } else if (gcount2 > 1) {
		call sprintf (Memc[errmsg], SZ_LINE, nogroup)
		call pargstr (name)
		call error (1, Memc[errmsg])
	    }
	}

	if (gcount2 == 0) {
	    if (ndim == 0 || ndim == dim1) {
		dim2 = dim1
		gcount2 = gcount1
	    } else if (ndim < dim1) {
		dim2 = ndim
		gcount2 = gcount1
		do idim = dim2+1, dim1 {
		    gcount2 = gcount2 * Meml[size+idim-1]
		    Meml[size+idim-1] = 1
		}
	    } else {
		dim2 = min (ndim, IM_MAXDIM)
		gcount2 = 1
		Meml[size+dim1] = gcount1
	    }

	} else if (gcount1 == gcount2) {
	    dim2 = dim1

	} else if (gcount2 == 1) {
	    if (ndim == 0) {
		dim2 = dim1 + 1
		Meml[size+dim1] = gcount1
	    } else if (ndim > dim1) {
		dim2 = min (ndim, IM_MAXDIM)
		Meml[size+dim1] = gcount1
	    } else {
		call sprintf (Memc[errmsg], SZ_LINE, nomatch)
		call pargstr (name)
		call error (1, Memc[errmsg])
	    }

	} else {
	    call sprintf (Memc[errmsg], SZ_LINE, nomatch)
	    call pargstr (name)
	    call error (1, Memc[errmsg])
	}

	# Rebuild image name from parts

	if (gcount2 == 1) {
	    call strcpy (IMG_NAME(img), Memc[fullname], SZ_FNAME)
	    call strcat (IMG_SECTION(img), Memc[fullname], SZ_FNAME)

	} else {
	    call sprintf (Memc[fullname], SZ_FNAME, "%s[%d/%d]%s")
	    call pargstr (IMG_NAME(img))
	    call pargi (gindex)
	    call pargi (gcount2)
	    call pargstr (IMG_SECTION(img))
	}

	# Check to see if new image is a pixel list

	if (chkexten (name, "pl")) {
	    IMG_PIXLIST(img) = YES
	    im2 = immap (Memc[fullname], NEW_FILE, NULL)
	    call pl_ssize (IM_PL(im2), dim2, Meml[size], 0)
	    call imaddi (im2, TYP_KWORD, type)

	} else {
	    IMG_PIXLIST(img) = NO
	    im2 = immap (Memc[fullname], NEW_FILE, NULL)
	    IM_PIXTYPE(im2) = type
	}

	# Set the pixel type and dimensions of the new image

	IM_NDIM(im2) = dim2
	IM_NPHYSDIM(im2) = dim2

	do idim = 1, IM_MAXDIM {
	    IM_LEN(im2,idim) = Meml[size+idim-1]
	    IM_SVLEN(im2,idim) = Meml[size+idim-1]
	}

	# Fill in fields of image descriptor

	IMG_DESCRIP(img) = im2
	IMG_MODE(img) = NEW_FILE
	IMG_GCOUNT(img) = gcount2

	call sfree (sp)
end

# NXT_IMAGE -- Open the next group in an image

int procedure nxt_image (img, gindex, line)

pointer img		# i: image descriptor
int	gindex		# u: index to current group
long	line[IM_MAXDIM]	# u: line to write
#--
pointer	sp, fullname, im, oldim

pointer	immap()

begin
	# See if there are any more groups

	gindex = gindex + 1
	if (gindex > IMG_GCOUNT(img))
	    return (NO)

	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (fullname, SZ_FNAME, TY_CHAR)

	# Unmap old image

	im = IMG_DESCRIP(img)
	if (IMG_MODE(img) == NEW_COPY) {
	    oldim = IM_OHDR(im)
	} else {
	    oldim = NULL
	}

	call imunmap (im)

	# If more groups, open the next group

	call sprintf (Memc[fullname], SZ_FNAME, "%s[%d]%s")
	call pargstr (IMG_NAME(img))
	call pargi (gindex)
	call pargstr (IMG_SECTION(img))

	IMG_DESCRIP(img) = immap (Memc[fullname], IMG_MODE(img), oldim)

	call amovkl (long(1), line, IM_MAXDIM)

	call sfree (sp)
	return (YES)
end


# OLD_IMAGE -- Open an existing image file

procedure old_image (name, mode, gindex, img)

char	name[ARB]	# i: name of image file
int	mode		# i: image access mode
int	gindex		# o: index to current group
pointer	img		# o: image descriptor
#--
int	gcount
pointer	sp, fullname, im

bool	chkexten()
int	imgeti(), btoi()
pointer	immap()

begin
	call smark (sp)
	call salloc (fullname, SZ_FNAME, TY_CHAR)

	# Allocate image structure

	call malloc (img, LEN_IMGSTRUCT, TY_INT)
	call malloc (IMG_NAMPTR(img), SZ_FNAME, TY_CHAR)
	call malloc (IMG_SECTPTR(img), SZ_FNAME, TY_CHAR)

	# Get information from image fields

	call parse_image (name, gindex, gcount, 
			 IMG_NAME(img), IMG_SECTION(img), SZ_FNAME)

	# Re-create image name and open it

	if (gcount == 0) {
	    # Don't use group number if not found in name
	    call strcpy (IMG_NAME(img), Memc[fullname], SZ_FNAME)
	    call strcat (IMG_SECTION(img), Memc[fullname], SZ_FNAME)

	} else {
	    call sprintf (Memc[fullname], SZ_FNAME, "%s[%d]%s")
	    call pargstr (IMG_NAME(img))
	    call pargi (gindex)
	    call pargstr (IMG_SECTION(img))
	}

	im = immap (Memc[fullname], mode, NULL)

	# Get group count from image 
	# if it was not specified in the name

	if (gcount == 0) {
	    iferr (gcount = imgeti (im, "GCOUNT")) {
		gcount = 1
	    } else if (IS_INDEFI (gcount)) {
		gcount = 1
	    }
	}

	# Fill in fields of image descriptor

	IMG_DESCRIP(img) = im
	IMG_PIXLIST(img) = btoi (chkexten (name, "pl"))
	IMG_MODE(img) = mode
	IMG_GCOUNT(img) = gcount

	call sfree (sp)
end

# PAR_IMAGE -- Get the parameters associated with an image

procedure par_image (img, type, nsize, gcount, size)

pointer img		# i: image descriptor
int	type		# o: image data type
int	nsize		# o: dimensionality of mask
int	gcount		# o: number of groups in mask
long	size[IM_MAXDIM]	# o: dimensions of mask
#--
int	idim
pointer	im

int	imgeti()

begin
	im = IMG_DESCRIP(img)

	if (IMG_PIXLIST(img) == NO) {
	    type = IM_PIXTYPE(im)
	} else {
	    iferr (type = imgeti (im, TYP_KWORD))
		type = IM_PIXTYPE(im)
	}

	nsize = IM_NDIM(im)
	gcount = IMG_GCOUNT(img)

	do idim = 1, IM_MAXDIM
	    size[idim] = IM_LEN(im,idim)

end

# PARSE_IMAGE -- Parse a image name into its component parts

procedure parse_image (name, gindex, gcount, root, section, maxch)

char	name[ARB]	# i: image name
int	gindex		# o: group index
int	gcount		# o: group count
char	root[ARB]	# o: image root
char	section[ARB]	# o: image section
int	maxch		# i: length of image section
#--
bool	star, valid
int	ifield, inum, num[2]
pointer	sp, ch, fields, field[NFIELD]

begin
	# Allocate memory to hold bracket delimeted fields in image name

	call smark (sp)
	call salloc (fields, (SZ_FNAME+1)*NFIELD, TY_CHAR)

	do ifield = 1, NFIELD
	    field[ifield] = fields + (ifield - 1) * (SZ_FNAME + 1)

	# Break image name into bracket delimeted fields

	call tp_break (name, Memc[fields], NFIELD, SZ_FNAME)
	call strcpy (Memc[field[1]], root, maxch)

	# Parse the second field as the group specifier

	inum = 0
	num[1] = 0
	num[2] = 0
	star = false
	valid = true

	for (ch = field[2]; Memc[ch] != EOS; ch = ch + 1) {
	    switch (Memc[ch]) {
	    case ' ':
		;
	    case '[':
		inum = 1
	    case ']':
		break
	    case '*':
		star = true
	    case '/':
		inum = inum + 1
	    default:
		if (IS_DIGIT(Memc[ch]) && ! star && (inum == 1 || inum == 2)) {
		    num[inum] = 10 * num[inum] + TO_INTEG(Memc[ch])

		} else {
		    star = false
		    valid = false
		    break
		}
	    }
	}

	# Extract the group index and count

	if (! valid) {
	    call strcpy (Memc[field[2]], section, maxch)
	    gindex = 1
	    gcount = 0

	} else {
	    section[1] = EOS

	    switch (inum) {
	    case 0:
		gindex = 1
		gcount = 0
	    case 1:
		gindex = max (1, num[1])
		gcount = 1
	    case 2:
		gindex = max (1, num[1])
		gcount = max (1, num[2])
	    }

	    if (star)
		gcount = 0

	}

	do ifield = 3, NFIELD
	    call strcat (Memc[field[ifield]], section, maxch)

end

# RD_IMAGE -- Read the next line from the image file

pointer procedure rd_image (type, img, gindex, line)

int	type		# i: data type of line
pointer img		# i: image descriptor
int	gindex		# u: index to current group
long	line[IM_MAXDIM]	# u: line to read
#--
int	status
pointer	vector

int	imgnld(), imgnli(), imgnlr(), nxt_image()

begin
	repeat {
	    # Read the next line from the instrument

	    switch (type) {
	    case TY_SHORT, TY_INT, TY_LONG, TY_USHORT:
		status = imgnli (IMG_DESCRIP(img), vector, line)
	    case TY_REAL:
		status = imgnlr (IMG_DESCRIP(img), vector, line)
	    case TY_DOUBLE:
		status = imgnld (IMG_DESCRIP(img), vector, line)
	    default:
		status = imgnld (IMG_DESCRIP(img), vector, line)
	    }

	    # If at EOF, get next group

	    if (status == EOF) {
		if (nxt_image (img, gindex, line) == NO) {
		    vector = NULL
		    status = OK
		}
	    }

	} until (status != EOF)

	return (vector)
end

# WRT_IMAGE -- Read a line to the image file

pointer procedure wrt_image (type, img, gindex, line)

int	type		# i: type of line
pointer img		# i: image descriptor
int	gindex		# u: index to current group
long	line[IM_MAXDIM]	# u: line to write
#--
int	status
pointer	vector

int	impnld(), impnli(), impnlr(), nxt_image()

begin
	repeat {
	    # Write the next line to the instrument

	    switch (type) {
	    case TY_SHORT, TY_INT, TY_LONG, TY_USHORT:
		status = impnli (IMG_DESCRIP(img), vector, line)
	    case TY_REAL:
		status = impnlr (IMG_DESCRIP(img), vector, line)
	    case TY_DOUBLE:
		status = impnld (IMG_DESCRIP(img), vector, line)
	    default:
		status = impnld (IMG_DESCRIP(img), vector, line)
	    }

	    # If at EOF, get next group

	    if (status == EOF) {
		if (nxt_image (img, gindex, line) == NO) {
		    vector = NULL
		    status = OK
		}
	    }

	} until (status != EOF)

	return (vector)
end


include <imhdr.h>

# These routines perform i/o to multigroup images, bad pixel lists,
# and tables, disguising the differences between them so that higher 
# level routines need not worry about the differences

# Originally written by Bernie Simon?
# Phil Hodge, 17-Oct-2002  Modify new_mask() to rely on filename extensions
#			or type of input mask, rather than using is_image().

define	LEN_MSKSTRUCT	(IM_MAXDIM+3)

define	MSK_DESCRIP	Meml[$1]		# image or table descriptor
define	MSK_IMAGE	Meml[$1+1]		# is mask an image?
define	MSK_GINDEX	Meml[$1+2]		# index to current group
define	MSK_LINE	Meml[$1+3]		# index to current line

# CLS_MASK -- Close mask and free associated structure

procedure cls_mask (mask)

pointer	mask		# i: mask descriptor
#--

begin
	# Free image or table descriptor

	if (MSK_IMAGE(mask) == YES)
	    call cls_image (MSK_DESCRIP(mask))
	else
	    call cls_table (MSK_DESCRIP(mask))

	# Free mask descriptor

	call mfree (mask, TY_LONG)
end

# DEL_MASK -- Delete mask and free associated structure

procedure del_mask (mask)

pointer	mask		# i: mask descriptor
#--

begin
	# Delete image or table

	if (MSK_IMAGE(mask) == YES)
	    call del_image (MSK_DESCRIP(mask))
	else
	    call del_table (MSK_DESCRIP(mask))

	# Free mask descriptor

	call mfree (mask, TY_LONG)
end

# OLD_MASK -- Open an existing mask file

pointer procedure old_mask (name, mode)

char	name[ARB]	# i: name of mask file
int	mode		# i: image access mode
#--
pointer	sp, errmsg, mask

string	badmode   "old_mask: cannot open old mask NEW_IMAGE or NEW_COPY"
string	notype    "Cannot determine mask file type (%s)"

int	is_image()

begin
	call smark (sp)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	# Can't open existing image with these modes

	if (mode == NEW_COPY || mode == NEW_FILE)
	    call error (1, badmode)

	# Allocate mask structure

	call malloc (mask, LEN_MSKSTRUCT, TY_LONG)
	call amovkl (long(1), MSK_LINE(mask), IM_MAXDIM)

	# Determine if mask is image or table

	MSK_IMAGE(mask) = is_image (name)

	# Open it accordingly

	if (MSK_IMAGE(mask) == YES) {
	    call old_image (name, mode, MSK_GINDEX(mask), MSK_DESCRIP(mask))

	} else if (MSK_IMAGE(mask) == NO) {
	    MSK_GINDEX(mask) = 1
	    call old_table (name, mode, MSK_DESCRIP(mask))

	} else {
	    call sprintf (Memc[errmsg], SZ_LINE, notype)
	    call pargstr (name)
	    call error (1, Memc[errmsg])
	}

	call sfree (sp)
	return (mask)
end

# NEW_MASK -- Open a new mask file

pointer procedure new_mask (name, ndim, omask)

char	name[ARB]	# i: name of mask file
int	ndim		# i: dimesionality of new mask (0 means same as old)
pointer	omask		# i: old mask descriptor
#--
pointer	sp, size, errmsg, mask

bool	chkexten()

string	notype    "Cannot determine mask file type (%s)"

begin
	call smark (sp)
	call salloc (errmsg, SZ_LINE, TY_CHAR)
	call salloc (size, IM_MAXDIM, TY_LONG)

	# Allocate mask structure

	call malloc (mask, LEN_MSKSTRUCT, TY_LONG)
	call amovkl (long(1), MSK_LINE(mask), IM_MAXDIM)

	# Determine if mask is image or table, and open it accordingly.
	# If we can't tell the type from the filename extension, use the
	# type of the old mask.

	if (chkexten (name, "??h") || chkexten (name, "pl"))
	    MSK_IMAGE(mask) = YES
	else if (chkexten (name, "tab"))
	    MSK_IMAGE(mask) = NO
	else
	    MSK_IMAGE(mask) = MSK_IMAGE(omask)

	if (MSK_IMAGE(mask) == YES) {
	    if (MSK_IMAGE(omask) == YES) {
		call cpy_image (name, ndim, MSK_DESCRIP(omask), 
				MSK_GINDEX(mask), MSK_DESCRIP(mask))

	    } else {
		call new_image (name, ndim, MSK_DESCRIP(omask),
				MSK_GINDEX(mask), MSK_DESCRIP(mask))
	    }

	} else if (MSK_IMAGE(mask) == NO){
	    if (MSK_IMAGE(omask) == YES) {
		MSK_GINDEX(mask) = 1
		call new_table (name, ndim, MSK_DESCRIP(omask), 
				MSK_DESCRIP(mask))
	    } else {
		MSK_GINDEX(mask) = 1
		call cpy_table (name, ndim, MSK_DESCRIP(omask), 
				MSK_DESCRIP(mask))
	    }

	} else {
	    call sprintf (Memc[errmsg], SZ_LINE, notype)
	    call pargstr (name)
	    call error (1, Memc[errmsg])
	}

	call sfree (sp)
	return (mask)
end

# RD_MASK -- Read the next line from the mask file

pointer procedure rd_mask (mask)

pointer mask		# i: mask descriptor
#--
pointer	vector

pointer	rd_image(), rd_table()

begin
	if (MSK_IMAGE(mask) == YES) {
	    vector = rd_image (TY_INT, MSK_DESCRIP(mask), 
			       MSK_GINDEX(mask), MSK_LINE(mask))
	} else {
	    vector = rd_table (MSK_DESCRIP(mask), MSK_GINDEX(mask),
			       MSK_LINE(mask))
	}

	return (vector)
end

# SZ_MASK -- Get the line length and number of lines in a mask file

procedure sz_mask (mask, length, nline)

pointer mask		# i: mask descriptor
int	length		# o: line length
int	nline		# o: number of lines
#--
int	type, nsize, gcount, idim
pointer	sp, size

begin
	call smark (sp)
	call salloc (size, IM_MAXDIM, TY_LONG)

	if (MSK_IMAGE(mask) == YES) {
	    call par_image (MSK_DESCRIP(mask), type, nsize, gcount, Meml[size])
	} else {
	    call par_table (MSK_DESCRIP(mask), type, nsize, gcount, Meml[size])
	}

	nline = gcount
	length = Meml[size]
	do idim = 2, nsize
	    nline = nline * Meml[size+idim-1]

	call sfree (sp)
end

# WRT_MASK -- Read the next line from the mask file

pointer procedure wrt_mask (mask)

pointer mask		# i: mask descriptor
#--
pointer	vector

pointer	wrt_image(), wrt_table()

begin
	if (MSK_IMAGE(mask) == YES) {
	    vector = wrt_image (TY_INT, MSK_DESCRIP(mask), 
				MSK_GINDEX(mask), MSK_LINE(mask))
	} else {
	    vector = wrt_table (MSK_DESCRIP(mask), MSK_GINDEX(mask),
				MSK_LINE(mask))
	}

	return (vector)
end


include <imhdr.h>

define	LEN_FIMSTRUCT	(IM_MAXDIM+3)

define	FIM_DESCRIP	Meml[$1]		# image descriptor
define	FIM_TYPE	Meml[$1+1]		# image type
define	FIM_GINDEX	Meml[$1+2]		# index to current group
define	FIM_LINE	Meml[$1+3]		# index to current line

define	FIM_LINPTR	($1+3)			# pointer to index to line

# BUP_FULLIMG -- Back up the line number in the image by one

procedure bup_fullimg (img)

pointer	img		# i: image descriptor
#--
int	type, nsize,gcount, idim
pointer	sp, size, dim

string	nobackup "bup_fullimg: cannot back up line number"

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (size, IM_MAXDIM, TY_LONG)

	# Get dimensions of the image

	call par_image (FIM_DESCRIP(img), type, nsize, gcount, Meml[size])

	# Back up line pointer

	dim = FIM_LINPTR(img)
	for (idim = 2; idim <= IM_MAXDIM; idim = idim + 1) {
	    dim = dim + 1

	    if (Meml[dim] == 1) {
		Meml[dim] = Meml[size+idim-1]
	    } else {
		Meml[dim] = Meml[dim] - 1
		break
	    }
	}

	# Check to see if number of dimensions was exceeded

	if (idim > IM_MAXDIM)
	    call error (1, nobackup)

	call sfree (sp)
end

# CLS_FULLIMG -- Close image and free associated structure

procedure cls_fullimg (img)

pointer	img		# i: image descriptor
#--

begin
	call cls_image (FIM_DESCRIP(img))
	call mfree (img, TY_LONG)
end

# OPN_FULLIMG -- Open an existing full image file

pointer procedure opn_fullimg (name, mode)

char	name[ARB]	# i: name of image file
int	mode		# i: image access mode
#--
int	nsize, gcount
pointer	img

string	badmode  "opn_fullimg: cannot open full image NEW_IMAGE or NEW_COPY"

begin
	# Can't open existing image with these modes

	if (mode == NEW_COPY || mode == NEW_FILE)
	    call error (1, badmode)

	# Allocate image structure

	call malloc (img, LEN_FIMSTRUCT, TY_LONG)

	# Fill in structure fields

	call old_image (name, mode, FIM_GINDEX(img), FIM_DESCRIP(img))

	call par_image (FIM_DESCRIP(img), FIM_TYPE(img), 
			nsize, gcount, FIM_LINE(img))

	call amovkl (long(1), FIM_LINE(img), IM_MAXDIM)

	return (img)
end

# RD_FULLIMG -- Read the next line from the image file

pointer procedure rd_fullimg (img)

pointer img		# i: image descriptor
#--
pointer	vector

pointer	rd_image()

begin
	vector = rd_image (FIM_TYPE(img), FIM_DESCRIP(img), 
			   FIM_GINDEX(img), FIM_LINE(img))
	return (vector)
end

# SZ_FULLIMG -- Get the line length and number of lines in a full image

procedure sz_fullimg (img, length, nline)

pointer img		# i: image descriptor
int	length		# o: line length
int	nline		# o: number of lines
#--
int	type, nsize, gcount, idim
pointer	sp, size

begin
	call smark (sp)
	call salloc (size, IM_MAXDIM, TY_LONG)

	call par_image (FIM_DESCRIP(img), type, nsize, gcount, Meml[size])

	nline = gcount
	length = Meml[size]
	do idim = 2, nsize
	    nline = nline * Meml[size+idim-1]

	call sfree (sp)
end

# TYP_FULLING -- Get the type of a full image file

procedure typ_fullimg (img, type)

pointer img		# i: image descriptor
int	type		# o: image type
#--

begin
	type = FIM_TYPE(img)
end

# WRT_FULLIMG -- Read the next line from the image file

pointer procedure wrt_fullimg (img)

pointer img		# i: image descriptor
#--
pointer	vector

pointer	wrt_image()

begin
	vector = wrt_image (FIM_TYPE(img), FIM_DESCRIP(img), 
			    FIM_GINDEX(img), FIM_LINE(img))
	return (vector)
end


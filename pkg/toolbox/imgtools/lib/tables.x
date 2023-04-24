include <imhdr.h>
include <tbset.h>
include <ctype.h>

define	LEN_TBLSTRUCT	(6+IM_MAXDIM)

define	TBL_DESCRIP	Memi[$1]		# table descriptor
define	TBL_BUFFER	Memi[$1+1]		# line buffer
define	TBL_INDEX	Memi[$1+2]		# current row number
define	TBL_DIRTY	Memi[$1+3]		# line buffer not flushed
define	TBL_DONE	Memi[$1+4]		# all lines in mask read
define	TBL_GINDEX	Memi[$1+5]		# group index of buffer
define	TBL_LINE	Memi[$1+5+$2]		# line index of buffer

define	GRP_CNAME	"GROUP"			# column name of group index
define	POS_CNAME	"PIX1"			# column names of position idx
define	LEN_CNAME	"LENGTH"		# column name of length info
define	VAL_CNAME	"VALUE"			# column name of value info

define	TYP_KWORD	"PIXTYPE"		# data type of mask
define	GRP_KWORD	"GROUPNUM"		# number of groups in mask
define	DIM_KWORD	"NAXES"			# number of axes in mask
define	LEN_KWORD	"AXLEN1"		# length of each axis

# CLS_TABLE -- Close table and free associated structure

procedure cls_table (tbl)

pointer	tbl		# i: table descriptor
#--
pointer	tp

errchk	imunmap

begin
	tp = TBL_DESCRIP(tbl)

	# Flush line buffer

	if (TBL_DIRTY(tbl) == YES)
	    call upd_table (tbl)

	# Free table descriptor

	call mfree (TBL_BUFFER(tbl), TY_INT)
	call mfree (tbl, TY_INT)

	# Close table

	if (tp != NULL)
	    call tbtclo (tp)

end

# CMP_TABLE -- Compare group index and line array to table row

int procedure cmp_table (tbl, gindex, line)

pointer	tbl		# i: Table descriptor
int	gindex		# i: Current group index
long	line[IM_MAXDIM]	# i: Current line
#--
int	row, len, value, idim, ndim, order
pointer	sp, tp, cp, cname

int	gstrcpy()

begin
	call smark (sp)
	call salloc (cname, SZ_COLNAME, TY_CHAR)

	tp = TBL_DESCRIP(tbl)
	row = TBL_INDEX(tbl)

	# Compare group indices. The result is the difference between
	# the value in the table and the input group index

	call tbcfnd (tp, GRP_CNAME, cp, 1)
	if (cp == NULL) {
	    order = 0
	} else {
	     call tbegti (tp, cp, row, value)
	     order = value - gindex
	}

	# Compare line number to table column from most to least
	# significant

	call hdr_table (tp, DIM_KWORD, ndim)
	len = gstrcpy (POS_CNAME, Memc[cname], SZ_COLNAME)

	for (idim = ndim; idim > 1 && order == 0; idim = idim - 1) {

	    Memc[cname+len-1] = TO_DIGIT(idim)
	    call colptr_table (tp, Memc[cname], YES, cp)
	    call tbegti (tp, cp, row, value)

	    order = value - line[idim]
	}

	call sfree (sp)
	return (order)
end

# COLPTR_TABLE -- Get a column pointer to a table

procedure colptr_table (tp, colname, crash, ptr)

pointer	tp		# i: table pointer
char	colname[ARB]	# i: column name
int	crash		# i: stop with error message if not found
pointer	ptr		# o: column pointer
#--
pointer	sp, tabname, errmsg

string	nocolumn  "Column %s not found in mask table (%s)"

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (tabname, SZ_FNAME, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	# Get column pointer from name
	# If not found, print error message

	call tbcfnd (tp, colname, ptr, 1)
	if (crash == YES && ptr == NULL) {
	    call tbtnam (tp, Memc[tabname], SZ_FNAME)

	    call sprintf (Memc[errmsg], SZ_LINE, nocolumn)
	    call pargstr (colname)
	    call pargstr (Memc[tabname])

	    call error (1, Memc[errmsg])
	}

	call sfree (sp)
end

# CPY_TABLE -- Create a new table which is a copy of an old one

procedure cpy_table (name, ndim, oldtbl, newtbl)

char	name[ARB]	# i: name of new table
int	ndim		# i: dimesionality of new mask (0 means same as old)
pointer	oldtbl		# i: old descriptor
pointer	newtbl		# o: new descriptor
#--
int	idim, dim1, dim2, gcount1, gcount2, type
pointer	sp, size, tp1, tp2

pointer	tbtopn()

begin
	call smark (sp)
	call salloc (size, IM_MAXDIM, TY_LONG)

	tp1 = TBL_DESCRIP(oldtbl)
	call par_table (oldtbl, type, dim1, gcount1, Meml[size])

	# Set the dimensions of the new table

	if (ndim == 0) {
	    tp2 = tbtopn (name, NEW_COPY, tp1)
	    call tbtcre (tp2)
	    call tbhcal (tp1, tp2)

	} else if (ndim == dim1) {
	    tp2 = tbtopn (name, NEW_COPY, tp1)
	    call tbtcre (tp2)
	    call tbhcal (tp1, tp2)

	} else if (ndim < dim1) {
	    dim2 = ndim
	    gcount2 = gcount1
	    do idim = dim2+1, dim1 {
		gcount2 = gcount2 * Meml[size+idim-1]
		Meml[size+idim-1] = 1
	    }

	    call make_table (name, type, dim2, gcount2, Meml[size], tp2)

	} else {
	    dim2 = min (ndim, IM_MAXDIM)
	    gcount2 = 1
	    Meml[size+dim1] = gcount1

	    call make_table (name, type, dim2, gcount2, Meml[size], tp2)
	}

	# Fill in fields of table descriptor

	call malloc (newtbl, LEN_TBLSTRUCT, TY_INT)
	call malloc (TBL_BUFFER(newtbl), Meml[size], TY_INT)

	TBL_DESCRIP(newtbl) = tp2
	TBL_INDEX(newtbl) = 1
	TBL_DIRTY(newtbl) = NO
	TBL_DONE(newtbl) = NO
	TBL_GINDEX(newtbl) = 1

	call amovki (1, TBL_LINE(newtbl,1), IM_MAXDIM)

	call sfree (sp)
end

# DEL_TABLE -- Delete table and free associated structure

procedure del_table (tbl)

pointer	tbl		# i: table descriptor
#--
pointer	tp, sp, tabname

errchk	tbtclo, delete

begin
	call smark (sp)
	call salloc (tabname, SZ_FNAME, TY_CHAR)

	tp = TBL_DESCRIP(tbl)

	# Close table

	if (tp != NULL) {
	    call tbtnam (tp, Memc[tabname], SZ_FNAME)

	    call tbtclo (tp)
	    call delete (Memc[tabname])
	}

	# Free table descriptor

	call mfree (TBL_BUFFER(tbl), TY_INT)
	call mfree (tbl, TY_INT)

	call sfree (sp)
end

# HDR_TABLE -- Get a header keyword value from a table

procedure hdr_table (tp, keyword, value)

pointer	tp		# i: table pointer
char	keyword[ARB]	# i: keyword name
int	value		# o: keyword value
#--
pointer	sp, tabname, errmsg

string	notfound  "Keyword %s not found in mask table header (%s)"

int	tbhgti()

begin
	call smark (sp)
	call salloc (tabname, SZ_FNAME, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	iferr {
	    value = tbhgti (tp, keyword)
	} then {
	    call tbtnam (tp, Memc[tabname], SZ_FNAME)

	    call sprintf (Memc[errmsg], SZ_FNAME, notfound)
	    call pargstr (keyword)
	    call pargstr (Memc[tabname])

	    call error (1, Memc[errmsg])
	}

	call sfree (sp)
end

# INC_TABLE -- Increment the line and group indices for a table

procedure inc_table (tbl, gindex, line, done)

pointer	tbl		# i: table descriptor
int	gindex		# u: group index
long	line[IM_MAXDIM]	# u: line index
int	done		# o: last line in table reached
#--
int	len, idim, ndim, maxidx, gcount
pointer	sp, tp, keyword

int	gstrcpy()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (keyword, SZ_KEYWORD, TY_CHAR)

	tp = TBL_DESCRIP(tbl)

	# Get the number of axes in the mask

	call hdr_table (tp, DIM_KWORD, ndim)

	# Increment the line index array, 
	# checking it against the mask dimensions

	done = YES
	len = gstrcpy (LEN_KWORD, Memc[keyword], SZ_KEYWORD)

	do idim = 2, ndim {
	    Memc[keyword+len-1] = TO_DIGIT(idim)
	    call hdr_table (tp, Memc[keyword], maxidx)

	    if (line[idim] >= maxidx) {
		line[idim] = 1
	    } else {
		done = NO
		line[idim] = line[idim] + 1
		break
	    }
	}

	# Increment the group count, when necessary

	if (done == YES) {
	    call hdr_table (tp, GRP_KWORD, gcount)

	    if (gindex < gcount) {
		done = NO
		gindex = gindex + 1
	    }
	}

	call sfree (sp)
end

# MAKE_TABLE -- Make a new table with the specified dimensions

procedure make_table (name, type, ndim, gcount, size, tp)

char	name[ARB]	# i: table name
int	type		# i: data type of table
int	ndim		# i: number of dimensions in mask
int	gcount		# i: number of groups in mask
long	size[IM_MAXDIM]	# i: length of each dimension
pointer	tp		# o: table pointer
#--
int	len, idim, length
pointer	sp, cp, cname, hname

string	units  ""
string	fmt    "%5d"

int	gstrcpy()
pointer	tbtopn()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (cname, SZ_COLNAME, TY_CHAR)
	call salloc (hname, SZ_KEYWORD, TY_CHAR)

	# Open table as new file

	tp = tbtopn (name, NEW_FILE, NULL)
	call tbpset (tp, TBL_MAXPAR, ndim+5)

	# Create table columns

	if (gcount > 1)
	    call tbcdef (tp, cp, GRP_CNAME, units, fmt, TY_INT, 1, 1)

	len = gstrcpy (POS_CNAME, Memc[cname], SZ_COLNAME)

	do idim = 1, ndim {
	    Memc[cname+len-1] = TO_DIGIT(idim)
	    call tbcdef (tp, cp, Memc[cname], units, fmt, TY_INT, 1, 1)
	}

	call tbcdef (tp, cp, LEN_CNAME, units, fmt, TY_INT, 1, 1)
	call tbcdef (tp, cp, VAL_CNAME, units, "", TY_INT, 1, 1)

	# Create table

	call tbtcre (tp)

	# Add header keywords

	call tbhadi (tp, TYP_KWORD, type)
	call tbhadi (tp, GRP_KWORD, gcount)
	call tbhadi (tp, DIM_KWORD, ndim)

	len = gstrcpy (LEN_KWORD, Memc[hname], SZ_KEYWORD)

	do idim = 1, ndim {
	    length = size[idim]
	    Memc[hname+len-1] = TO_DIGIT(idim)
	    call tbhadi (tp, Memc[hname], length)
	}

	call sfree (sp)
end

# NEW_TABLE -- Create a new table which is a copy of an image

procedure new_table (name, ndim, img, tbl)

char	name[ARB]	# i: name of new table
int	ndim		# i: dimesionality of new mask (0 means same as old)
pointer	img		# i: image descriptor
pointer	tbl		# o: table descriptor
#--
int	idim, dim1, dim2, gcount1, gcount2, type
pointer	sp, size, tp

begin
	call smark (sp)
	call salloc (size, IM_MAXDIM, TY_LONG)

	# Set the dimensions of the new table

	call par_image (img, type, dim1, gcount1, Meml[size])

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

	call make_table (name, type, dim2, gcount2, Meml[size], tp)


	# Fill in fields of table descriptor

	call malloc (tbl, LEN_TBLSTRUCT, TY_INT)
	call malloc (TBL_BUFFER(tbl), Meml[size], TY_INT)

	TBL_DESCRIP(tbl) = tp
	TBL_INDEX(tbl) = 1
	TBL_DIRTY(tbl) = NO
	TBL_DONE(tbl) = NO
	TBL_GINDEX(tbl) = 1

	call amovki (1, TBL_LINE(tbl,1), IM_MAXDIM)

	call sfree (sp)
end

# OLD_TABLE -- Open an existing table file

procedure old_table (name, mode, tbl)

char	name[ARB]	# i: name of table file
int	mode		# i: table access mode
pointer	tbl		# o: table descriptor
#--
int	nbuf
pointer	tp

pointer	tbtopn()

begin
	tp = tbtopn (name, mode, NULL)
	call hdr_table (tp, LEN_KWORD, nbuf)

	# Fill in fields of table descriptor

	call malloc (tbl, LEN_TBLSTRUCT, TY_INT)
	call malloc (TBL_BUFFER(tbl), nbuf, TY_INT)

	TBL_DESCRIP(tbl) = tp
	TBL_INDEX(tbl) = 1
	TBL_DIRTY(tbl) = NO
	TBL_DONE(tbl) = NO
	TBL_GINDEX(tbl) = 1

	call amovki (1, TBL_LINE(tbl,1), IM_MAXDIM)

end

# PAR_TABLE -- Get the parameters associated with a table mask

procedure par_table (tbl, type, nsize, gcount, size)

pointer tbl		# i: table descriptor
int	type		# i: mask data type
int	nsize		# o: dimensionality of mask
int	gcount		# o: number of groups in mask
long	size[IM_MAXDIM]	# o: dimensions of mask
#--
int	len, idim, length
pointer	sp, tp, keyword

int	gstrcpy()

begin
	call smark (sp)
	call salloc (keyword, SZ_KEYWORD, TY_CHAR)

	tp = TBL_DESCRIP(tbl)

	call hdr_table (tp, TYP_KWORD, type)
	call hdr_table (tp, DIM_KWORD, nsize)
	call hdr_table (tp, GRP_KWORD, gcount)

	call amovki (1, size, IM_MAXDIM)
	len = gstrcpy (LEN_KWORD, Memc[keyword], SZ_KEYWORD)

	do idim = 1, nsize {
	    Memc[keyword+len-1] = TO_DIGIT(idim)
	    call hdr_table (tp, Memc[keyword], length)
	    size[idim] = length
	}

	call sfree (sp)
end

# RD_TABLE -- Read the next line from the table file

pointer procedure rd_table (tbl, gindex, line)

pointer tbl		# i: table descriptor
int	gindex		# u: current group index
long	line[IM_MAXDIM]	# u: line to read
#--
int	nrow, nbuf, offset, length, value
pointer	tp, dim, val, len, buf

int	tbpsta(), cmp_table()

begin
	# Don't do anything if at end of table

	if (TBL_DONE(tbl) == YES)
	    return (NULL)

	# Get mask dimensions from table

	tp = TBL_DESCRIP(tbl)
	nrow = tbpsta (tp, TBL_NROWS)
	call hdr_table (tp, LEN_KWORD, nbuf)

	call colptr_table (tp, POS_CNAME, YES, dim)
	call colptr_table (tp, LEN_CNAME, NO, len)
	call colptr_table (tp, VAL_CNAME, YES, val)

	# Search for the current line

	while (TBL_INDEX(tbl) <= nrow && cmp_table (tbl, gindex, line) < 0)
	    TBL_INDEX(tbl) = TBL_INDEX(tbl) + 1

	# Use information in table to fill buffer with mask values

	buf = TBL_BUFFER(tbl)
	call amovki (0, Memi[buf], nbuf)

	while (TBL_INDEX(tbl) <= nrow && cmp_table (tbl, gindex, line) == 0) {
	    call tbegti (tp, dim, TBL_INDEX(tbl), offset)
	    call tbegti (tp, val, TBL_INDEX(tbl), value)

	    if (len == NULL) {
		Memi[buf+offset-1] = value
	    } else {
		call tbegti (tp, len, TBL_INDEX(tbl), length)
		call amovki (value, Memi[buf+offset-1], length)
	    }

	    TBL_INDEX(tbl) = TBL_INDEX(tbl) + 1
	}

	# Increment line number and group index
	
	call inc_table (tbl, gindex, line, TBL_DONE(tbl))
	return (buf)

end

# UPD_TABLE -- Write the buffer contents to the mask table

procedure upd_table (tbl)

pointer tbl		# i: table descriptor
#--
int	nbuf, ibuf, jbuf, len, last, length, nrow, nval, idim, ndim
pointer	sp, tp, colptr, value, keyword, buf

int	gstrcpy(), tbpsta()

begin
	call smark (sp)
	call salloc (colptr, IM_MAXDIM+3, TY_INT)
	call salloc (value, IM_MAXDIM+3, TY_INT)
	call salloc (keyword, SZ_KEYWORD, TY_CHAR)

	TBL_DIRTY(tbl) = NO

	# Get buffer dimensions from table

	tp = TBL_DESCRIP(tbl)
	call hdr_table (tp, LEN_KWORD, nbuf)

	# Check for null row

	buf = TBL_BUFFER(tbl)
	for (jbuf = 1; jbuf <= nbuf; jbuf = jbuf + 1) {
	    if (Memi[buf+jbuf-1] != 0) {
		last = Memi[buf+jbuf-1]
		length = 1
		break
	    }
	}
	
	if (jbuf > nbuf) {
	    call sfree (sp)
	    return
	}

	nrow = tbpsta (tp, TBL_NROWS)
	call hdr_table (tp, DIM_KWORD, ndim)

	call colptr_table (tp, GRP_CNAME, NO, Memi[colptr])
	if (Memi[colptr] == NULL) {
	    nval = 0
	} else {
	    Memi[value] = TBL_GINDEX(tbl)
	    nval = 1
	}

	len = gstrcpy (POS_CNAME, Memc[keyword], SZ_KEYWORD)
	do idim = ndim, 1, -1 {
	    Memc[keyword+len-1] = TO_DIGIT(idim)

	    call colptr_table (tp, Memc[keyword], YES, Memi[colptr+nval])
	    Memi[value+nval] = TBL_LINE(tbl,idim)

	    nval = nval + 1
	}

	call colptr_table (tp, LEN_CNAME, YES, Memi[colptr+nval])
	nval = nval + 1

	call colptr_table (tp, VAL_CNAME, YES, Memi[colptr+nval])
	nval = nval + 1

	do ibuf = jbuf+1, nbuf {
	    if (Memi[buf+ibuf-1] == last) {
		if (last != 0)
		    length = length + 1

	    } else if (last == 0) {
		last = Memi[buf+ibuf-1]
		length = 1
		jbuf = ibuf

	    } else {
		nrow = nrow + 1

		Memi[value+nval-3] = jbuf
		Memi[value+nval-2] = length
		Memi[value+nval-1] = last
		call tbrpti (tp, Memi[colptr], Memi[value], nval, nrow)

		last = Memi[buf+ibuf-1]
		if (last == 0) {
		    length = 0
		} else {
		    length = 1
		    jbuf = ibuf
		}
	    }
	}

	if (last != 0) {
	    nrow = nrow + 1

	    Memi[value+nval-3] = jbuf
	    Memi[value+nval-2] = length
	    Memi[value+nval-1] = last
	    call tbrpti (tp, Memi[colptr], Memi[value], nval, nrow)
	}

	call sfree (sp)
end

# WRT_TABLE -- Write a line to the table file

pointer procedure wrt_table (tbl, gindex, line)

pointer tbl		# i: table descriptor
int	gindex		# u: current group index
long	line[IM_MAXDIM]	# u: line to write
#--
int	idim

begin
	# Flush the buffer if it has already been written to

	if (TBL_DIRTY(tbl) == YES)
	    call upd_table (tbl)

	# Return a null pointer if the last line in the table
	# has been reached. Otherwise, copy the group and line
	# indexes into the structure, then increment them, then
	# return the buffer address

	if (TBL_DONE(tbl) == YES) {
	    return (NULL)

	} else {
	    TBL_DIRTY(tbl) = YES
	    TBL_GINDEX(tbl) = gindex
	    do idim = 1, IM_MAXDIM
		TBL_LINE(tbl,idim) = line[idim]

	    call inc_table (tbl, gindex, line, TBL_DONE(tbl))
	    return (TBL_BUFFER(tbl))
	}

end


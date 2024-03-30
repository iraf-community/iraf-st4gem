include	<finfo.h>
include <tbset.h>
include "libsynphot.h"

define	LEN_GRFSTRUCT	10

define	GRF_VALID	Memi[$1]	# address of structure, for validation
define	GRF_NROWS	Memi[$1+1]	# length of structure arrays
define	GRF_MTIME	Memi[$1+2]	# time graph table was last modified
define  GRF_OLDSTYLE    Memi[$1+3]      # if true, there is no thermal component column
define	GRF_TABLE	Memi[$1+4]	# graph table name
define	GRF_COMPARRAY	Memi[$1+5]	# optical component ids
define	GRF_KEYARRAY	Memi[$1+6]	# observation keywords
define	GRF_INARRAY	Memi[$1+7]	# graph input nodes
define	GRF_OUTARRAY	Memi[$1+8]	# graph output nodes
define  GRF_THMLARRAY   Memi[$1+9]      # thermal component ids

# MAPGRAF -- Load the instrument graph table into memory

pointer procedure mapgraf (graphtab)

char	graphtab[ARB]	# i: Graph table name
#--
char	star
int	mtime, nrow, irow, jrow, krow, innode, outnode
long	fi[LEN_FINFO]
pointer	sp, graph, compid, thcompid, keywrd, oldgrf, table, grf
pointer	tp, cp, flag, index, array, ptr1, ptr2
pointer tempid


data	star	/ '*' /
data	oldgrf	/ NULL /

string	compcol   "COMPNAME"
string  thermcol  "THCOMPNAME"
string	keycol    "KEYWORD"
string	inodecol  "INNODE"
string	onodecol  "OUTNODE"
string	nofiles1  "* Synphot data files may not be installed *\n"
string	nofiles2  "* See Appendix C of Installation Manual   *\n"
string	nofiles3  "* for further information on installation *\n" 
string	stars     "*******************************************\n"
string	nullgraph "Graph table is missing input node"
string	nograph   "Cannot access graph table"

bool	streq()
int	finfo(), tbpsta(), stridx()
pointer	tbtopn()

int     ordgraf()
extern	ordgraf
errchk	synphoterr, syncolptr, tbtopn, tbegti, tbegtt

begin


	# Get most recent graph table matching pattern

	call smark (sp)
	call salloc (graph, SZ_FNAME, TY_CHAR)
	call salloc (compid, SZ_COMPID, TY_CHAR)
	call salloc (thcompid, SZ_COMPID, TY_CHAR)
	call salloc (keywrd, SZ_KEYWRD, TY_CHAR)

	call lastfile (graphtab, Memc[graph], SZ_FNAME)

	# Check to see whether graph table has already been opened

	if (oldgrf == NULL) {
	    mtime = 0

	} else if (GRF_VALID(oldgrf) != oldgrf) {
	    mtime = 0

	} else {
	    table = GRF_TABLE(oldgrf)

	    if (streq (Memc[table], Memc[graph])) {
		if (finfo (Memc[graph], fi) == ERR) {
		    if (stridx (star, Memc[graph]) > 0) {
			call eprintf (stars)
			call eprintf (nofiles1)
			call eprintf (nofiles2)
			call eprintf (nofiles3)
			call eprintf (stars)
		    }

		    call synphoterr (nograph, Memc[graph])
		}

		if (FI_MTIME(fi) != GRF_MTIME(oldgrf)) {
		    mtime = FI_MTIME(fi)
		} else {
		    call sfree (sp)
		    grf = oldgrf
#call eprintf("this is the shortcut - new table not read \n")
		    return (grf)
		}
	    }
	}

#call eprintf("new table is read \n")

	# Open graph table and get size

	iferr {
	    tp = tbtopn (Memc[graph], READ_ONLY, NULL)
	    nrow = tbpsta (tp, TBL_NROWS)
	} then {
	    if (stridx (star, Memc[graph]) > 0) {
		call eprintf (stars)
		call eprintf (nofiles1)
		call eprintf (nofiles2)
		call eprintf (nofiles3)
		call eprintf (stars)
	    }

	    call synphoterr (nograph, Memc[graph])
	}

	# Allocate memory for structure

	if (oldgrf != NULL)
	    call unmapgraf (oldgrf)

	call malloc (grf, LEN_GRFSTRUCT, TY_STRUCT)
	call malloc (GRF_TABLE(grf), SZ_FNAME, TY_CHAR)
	call malloc (GRF_OLDSTYLE(grf), SZ_INT, TY_INT)
	call malloc (GRF_COMPARRAY(grf), nrow*(SZ_COMPID+1), TY_CHAR)
	call malloc (GRF_THMLARRAY(grf), nrow*(SZ_COMPID+1), TY_CHAR)
	call malloc (GRF_KEYARRAY(grf), nrow*(SZ_KEYWRD+1), TY_CHAR)
	call malloc (GRF_INARRAY(grf), nrow, TY_INT)
	call malloc (GRF_OUTARRAY(grf), nrow, TY_INT)

	call salloc (flag, nrow, TY_BOOL)
	call salloc (index, nrow, TY_INT)

# see if it's an old style table
	Memi[GRF_OLDSTYLE(grf)] = NO
        call tbcfnd (tp, thermcol, tempid, 1)
#	call eprintf("therm Column id is %d \n")
#	call pargi(thcompid)
	if (tempid <= 0) {
            Memi[GRF_OLDSTYLE(grf)] = YES
            call eprintf("Warning, you are using an old style graph table.\n")
	    call eprintf("No thermal background calculations will be supported.\n \n")
        } 


	# Initialize the other scalars in the structure
#call eprintf("oldstyle var is %b \n")
#call pargb(Memi[GRF_OLDSTYLE(grf)])


	GRF_VALID(grf) = grf
	GRF_NROWS(grf) = nrow
	table = GRF_TABLE(grf)
	call strcpy (Memc[graph], Memc[table], SZ_FNAME)

	if (mtime != 0) {
	    GRF_MTIME(grf) = mtime
	} else if (finfo (Memc[graph], fi) == ERR) {
	    call synphoterr (nograph, Memc[graph])
	} else {
	    GRF_MTIME(grf) = FI_MTIME(fi)
	}

	# Read in the graph table

#call eprintf("pre read \n")
	array = GRF_COMPARRAY(grf)
	call syncolptr (tp, compcol, 1, cp)
	call tbcgtt (tp, cp, Memc[array], Memb[flag], SZ_COMPID, 1, nrow)

	if (Memi[GRF_OLDSTYLE(grf)] == NO) {
	array = GRF_THMLARRAY(grf)
	call syncolptr (tp, thermcol, 5, cp)
	call tbcgtt (tp, cp, Memc[array], Memb[flag], SZ_COMPID, 1, nrow)
	} 

	array = GRF_KEYARRAY(grf)
	call syncolptr (tp, keycol, 2, cp)
	call tbcgtt (tp, cp, Memc[array], Memb[flag], SZ_KEYWRD, 1, nrow)

	array = GRF_OUTARRAY(grf)
	call syncolptr (tp, onodecol, 4, cp)
	call tbcgti (tp, cp, Memi[array], Memb[flag], 1, nrow)

	array = GRF_INARRAY(grf)
	call syncolptr (tp, inodecol, 3, cp)
	call tbcgti (tp, cp, Memi[array], Memb[flag], 1, nrow)


	do irow = 1, nrow {
	    if (Memb[flag])
		call synphoterr (nullgraph, Memc[graph])
	}

#call eprintf("pre sort \n")
	# Sort the graph table
	call iota (Memi[index], nrow)
	call gqsort (Memi[index], nrow, ordgraf, array)

	# Rearrange arrays in sorted order. The index array forms one or 
        # more cycles. The elements of these cycles are moved to the 
        # correct location and the index array is updated to indicate
        # which rows have been moved.

	do irow = 1, nrow {
	    if (Memi[index+irow-1] == irow)
		next

            # Move first element in cycle to temporary variables

	    ptr1 = GRF_COMPARRAY(grf) + (irow - 1) * (SZ_COMPID + 1)
	    call strcpy (Memc[ptr1], Memc[compid], SZ_COMPID)

	    if (Memi[GRF_OLDSTYLE(grf)] == NO) {
	    ptr1 = GRF_THMLARRAY(grf) + (irow - 1) * (SZ_COMPID + 1)
	    call strcpy (Memc[ptr1], Memc[thcompid], SZ_COMPID)
	    }

	    ptr1 = GRF_KEYARRAY(grf) + (irow - 1) * (SZ_KEYWRD + 1)
	    call strcpy (Memc[ptr1], Memc[keywrd], SZ_KEYWRD)

	    ptr1 = GRF_INARRAY(grf) + irow - 1
	    innode = Memi[ptr1]

	    ptr1 = GRF_OUTARRAY(grf) + irow - 1
	    outnode = Memi[ptr1]

            # Exchange elements in cycle

	    jrow = irow
	    while (Memi[index+jrow-1] != irow) {
		krow = Memi[index+jrow-1]

		ptr1 = GRF_COMPARRAY(grf) + (krow - 1) * (SZ_COMPID + 1)
		ptr2 = GRF_COMPARRAY(grf) + (jrow - 1) * (SZ_COMPID + 1)
		call strcpy (Memc[ptr1], Memc[ptr2], SZ_COMPID)

		if (Memi[GRF_OLDSTYLE(grf)] == NO) {
		ptr1 = GRF_THMLARRAY(grf) + (krow - 1) * (SZ_COMPID + 1)
		ptr2 = GRF_THMLARRAY(grf) + (jrow - 1) * (SZ_COMPID + 1)
		call strcpy (Memc[ptr1], Memc[ptr2], SZ_COMPID)
		}

		ptr1 = GRF_KEYARRAY(grf) + (krow - 1) * (SZ_KEYWRD + 1)
		ptr2 = GRF_KEYARRAY(grf) + (jrow - 1) * (SZ_KEYWRD + 1)
		call strcpy (Memc[ptr1], Memc[ptr2], SZ_KEYWRD)

		ptr1 = GRF_INARRAY(grf) + krow - 1
		ptr2 = GRF_INARRAY(grf) + jrow - 1
		Memi[ptr2] = Memi[ptr1]

		ptr1 = GRF_OUTARRAY(grf) + krow - 1
		ptr2 = GRF_OUTARRAY(grf) + jrow - 1
		Memi[ptr2] = Memi[ptr1]

		Memi[index+jrow-1] = jrow
		jrow = krow
	    }

#call eprintf("mid sort\n")

            # Move temporary variable into last element in cycle

	    ptr1 = GRF_COMPARRAY(grf) + (jrow - 1) * (SZ_COMPID + 1)
	    call strcpy (Memc[compid], Memc[ptr1], SZ_COMPID)

	    if (Memi[GRF_OLDSTYLE(grf)] == NO) {
	    ptr1 = GRF_THMLARRAY(grf) + (jrow - 1) * (SZ_COMPID + 1)
	    call strcpy (Memc[thcompid], Memc[ptr1], SZ_COMPID)
	    }

	    ptr1 = GRF_KEYARRAY(grf) + (jrow - 1) * (SZ_KEYWRD + 1)
	    call strcpy (Memc[keywrd], Memc[ptr1], SZ_KEYWRD)

	    ptr1 = GRF_INARRAY(grf) + jrow - 1
	    Memi[ptr1] = innode 

	    ptr1 = GRF_OUTARRAY(grf) + jrow - 1
	    Memi[ptr1] = outnode

	    Memi[index+jrow-1] = jrow
	}


	# Free temporary arrays and close table

	call tbtclo (tp)
	call sfree (sp)

	oldgrf = grf
	return (grf)
end

# ORDGRAF -- Return order of two elements in array

int procedure ordgraf (array, i1, i2)

pointer	array		# i: pointer to array
int	i1		# i: index to first element
int	i2		# i: index to second element
#--

begin
	return (Memi[array+i1-1] - Memi[array+i2-1])
end

# RDGRAFI -- Read integer value from graph structure

procedure rdgrafi (grf, icol, irow, value)

pointer	grf		# i: pointer to graph structure
int	icol		# i: column to read
int	irow		# i: row to read
int	value		# o: output value
#--
pointer	ptr

string	badcolumn  "Illegal column number in graph table"

begin
	switch (icol) {
	case GRF_INNODE:
	    ptr = GRF_INARRAY(grf) + irow - 1
	case GRF_OUTNODE:
	    ptr = GRF_OUTARRAY(grf) + irow - 1
	default:
	    call synphoterr (badcolumn, "rdgrafi")
	}

	value = Memi[ptr]
end

# RDGRAFT -- Read text value from graph structure

procedure rdgraft (grf, icol, irow, str, maxch)

pointer	grf		# i: pointer to graph structure
int	icol		# i: column to read
int	irow		# i: row to read
char	str[ARB]	# o: output string
int	maxch		# i:maximum length of output string
#--
pointer	ptr

string	badcolumn  "Illegal column number in graph table"
string  oldstyle   " has no thermal component column"
char    osmsg[SZ_FNAME+35]
begin

	switch (icol) {
	case GRF_COMPID:
	    ptr = GRF_COMPARRAY(grf) + (irow - 1) * (SZ_COMPID + 1)
	case GRF_THMLID:
	    if (Memi[GRF_OLDSTYLE(grf)] == NO) {
	      ptr = GRF_THMLARRAY(grf) + (irow - 1) * (SZ_COMPID + 1)
	    } 
            else {
	      call strcpy(Memc[GRF_TABLE(grf)], osmsg, SZ_FNAME)
	      call strcat(oldstyle, osmsg, SZ_FNAME+35]
	      call eprintf(" oldstyle var is %b \n")
	      call pargb(Memi[GRF_OLDSTYLE(grf)])
	      call synphoterr (osmsg, "rdgraft")
            }
	case GRF_KEYWRD:
	    ptr = GRF_KEYARRAY(grf) + (irow - 1) * (SZ_KEYWRD + 1)
	default:
	    call synphoterr (badcolumn, "rdgraft")
	}

	call strcpy (Memc[ptr], str, maxch)
end

# SIZGERAF -- Return number of rows in graph

int procedure sizegraf (grf)

pointer	grf		# u: pointer to graph structure
#--

begin
	return (GRF_NROWS(grf))
end

# UNMAPGRAF -- Free memory used by graph structure

procedure unmapgraf (grf)

pointer	grf		# u: pointer to graph structure
#--

begin
	call mfree (GRF_COMPARRAY(grf), TY_CHAR)
	call mfree (GRF_THMLARRAY(grf), TY_CHAR)
	call mfree (GRF_KEYARRAY(grf), TY_CHAR)
	call mfree (GRF_OUTARRAY(grf), TY_INT)
	call mfree (GRF_INARRAY(grf), TY_INT)
	call mfree (GRF_TABLE(grf), TY_CHAR)

	call mfree (grf, TY_STRUCT)
	grf = NULL
end


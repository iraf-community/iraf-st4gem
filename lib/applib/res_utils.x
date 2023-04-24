include	<pattern.h>
include	<ctype.h>
include	<time.h>
include	<tbset.h>
include	<foc.h>

define	SYNTAX		1
define	BOUNDS		2
define	PUTNULL		11
define	MAX_STACK	8
define	BLANK		' '
define	DELIM		','
define  ESCAPE		'\\'
define	NEGCHAR		'~'

# D. Giaretta, 01-June-1987	Original code
# D. Giaretta, 01-June-1988	Add error trapping to rs_open and xcoord_type
# Phil Hodge,  13-Mar-1989	Pass a char instead of an int to stridx.
# Phil Hodge,  19-Sep-1989	rs_open:  remove iferr;
#			rs_r_init:  check for column not found;
#			xcoord_type:  use nint, check for column not found,
#			and assign nrows,ncols = INDEFI.
# Phil Hodge,   3-Jan-1990	rs_readr:  change first test on ent_num from
#			(ent_num <= 0) to (ent_num < 0);
#			don't take last entry if no match.
# Phil Hodge, 20-Mar-1992	rs_readr:  check both X & Y for indef position.
# Phil Hodge, 13-Sep-1993	xcoord_type:  if coordfile is STDIN, read it
#			as a text file rather than a table; check this first.

.help res_utils
.nf _________________________________________________________________
RES_UTILS -- Routines to access reseau files, including routines
to initialise pointers and to read or write reseau entries into 
a reseau file

pointer rp

rp = rs_open( resfile, accmode, template)
		open reseau file, return pointer to file structure

call rs_close( rp )
		close reseau file

call rs_new( num_entries, ncols, rows, cols, datatype, rp)
		create new reseau file

int = rs_readr ( rp, respt, ent_name, ent_num )
		read reseau file - assumes rs_open called

call rs_write ( rp, respt, row_num)
		write reseau entry - assumes rs_open called

call rs_insert ( rpin, rp, respt )
		insert enrty into reseau file in correct alphanumeric order

call rs_cnvtime( respt, datetime, SZ_R_TIME)
		returns date/time string in datetime

call rs_alloc( rp, datatype, rentp ) 
		allocate space for entry information structure

bool 	= rs_check( rp, nrows, ncols)
		do two res grids have equal number of points?

int 	= rs_num_rows(rp)
		return the number of reseau entries in the file

call rs_t_update( respt, tinfo, len)
		adds tracking info into reseau entry info block

long 	= xcoord_type( coordfile, entry, xcol, ycol, xcoord, ycoord )
		given a file name, figure out if it is a text file or
		a table and then read the data into x and y arrays.
		returns number of coordinate pairs. If table file will 
		request the names "xcol" and "ycol" to read the x, y coords.
		Only non-INDEF coordinates are returned

.ih
Entry template package

This package contains subroutines to expand a entry name template into
an array of entry pointers which match the template.  The template is a
list of entry patterns separated by commas or whitespace.  The entry
pattern is either a entry name, a file name containing a list of entry
names, or a pattern using the usual IRAF pattern matching syntax. For
example, the string

	a[1-9], b, time*, @entry.lis

would be expanded as the entry names a1 through a9, b, any entry name
beginning with "time", and all the entry names in the file entry.lis.
If the entry template is entirely whitespace, the array of entry pointers
will include all the entrys in the table, as this seems the most reasonable
default. If the first non-white character of a template is the negation 
character (~), the array of entry pointers of entries matched by that
template will be marked as not to be used, although this may be over-ridden
by another template. The negation character only has this meaning as the 
first character in the entry template, and is interpreted as part of an 
entry name if found later in the template or in a file.

pointer = rs_template( rp, entryname)
		template expander for name "entryname"

int = rs_list(rlist)
		return entry number of next entry wanted, rlist returned from 
					rs_template

.endhelp ____________________________________________________________

# RS_OPEN -- open reseau file

pointer procedure rs_open( resfile, accmode, template)

char	resfile[SZ_FNAME]	# i: res. file name
int	accmode			# i: access mode
pointer	template		# i: template pointer
#--

pointer	tbtopn(), tp
pointer	rp
pointer	sp, nullc
int	num_tabrows
int	rows, cols
int	rs_num_rows(), numrows
errchk	tbtopn, tbrptt, tbrpti, tbtcre, tbhcal
errchk	rs_r_init

begin

	tp = NULL
	rp = NULL

	switch ( accmode) {
	case READ_ONLY, READ_WRITE, WRITE_ONLY:
	    if ( template != NULL)
		tp = tbtopn( resfile, accmode, RES_FILE(template) ) 
	    else
		tp = tbtopn( resfile, accmode, NULL ) 
	    call rs_r_init ( tp, rp)
	    numrows = rs_num_rows(rp)
	    # read all entry names into array
	    call salloc( RES_ALL_NAMES(rp), 
			 numrows*(SZ_R_ENTRY+1),TY_CHAR )
	    call smark(sp)
	    call salloc( nullc, num_tabrows, TY_BOOL)
	    call tbcgtt( tp, RES_COLPTR(rp, RES_P_ENTRY), 
				Memc[RES_ALL_NAMES(rp)], 
				Memb[nullc], SZ_R_ENTRY, 1, numrows)
	    call sfree ( sp )

	case NEW_FILE:
	    tp = tbtopn( resfile, accmode, NULL ) 
	    call tbpset( tp, TBL_MAXPAR, RES_USERPARS)
	    call salloc( rp, LEN_RES_INFO, TY_LONG)
	    RES_FILE(rp)  = tp

	case NEW_COPY:
	    tp = tbtopn( resfile, accmode, RES_FILE(template) ) 
	    # decide how many row and columns to make
	    rows = RES_ROWSPERENTRY*rs_num_rows(template)
	    call tbpset( tp, TBL_ALLROWS, rows )
	    call tbtcre(tp)

	    rows = RES_NROWS( template)
	    cols = RES_NCOLS( template)

	    # copy all header info parameters
	    call tbhcal( RES_FILE(template), tp)

	    call rs_r_init( tp, rp)
	    RES_FILE(rp)   = tp
	    RES_EXTENT(rp) = RES_EXTENT(template)
	    RES_NROWS( rp) = RES_NROWS( template)
	    RES_NCOLS( rp) = RES_NCOLS( template)

	default:
	    call error( 0, " illegal access mode for reseau file")
	}

	return ( rp )
end

# RS_R_INIT -- initialise column pointers to read/write reseau file

procedure rs_r_init ( tp, rp )

pointer	tp		# i: table pointer
pointer	rp		# o: pointer to reseau file structure
#--

pointer	colptr[NUM_R_INFOCOLS], tbcnum(), colpt
int	nrows, ncols, tbhgti()
char	i_colnames[SZ_COLNAME, NUM_R_INFOCOLS]
char	colname[SZ_COLNAME], colfmt[SZ_COLFMT], colunits[SZ_COLUNITS]
int	coldatatype, collendata, collenfmt
char	mess[SZ_RES_MESS]

int	i, strncmp(), ctoi(), ip, ival, colnum
int	k
int	tbpsta(), tabcols

errchk	tbhgti, tbcnum, tbpsta, tbcinf

begin

	call strcpy( "NAME",      i_colnames[1, RES_P_ENTRY], SZ_COLNAME)
	call strcpy( "TRACKING",  i_colnames[1, RES_P_TRACKING],  SZ_COLNAME)
	call strcpy( "DATE",      i_colnames[1, RES_P_DATE],  SZ_COLNAME)
	call strcpy( "NROWS",     i_colnames[1, RES_P_NROWS],  SZ_COLNAME)
	call strcpy( "NCOLS",     i_colnames[1, RES_P_NCOLS],  SZ_COLNAME)
	call strcpy( "SFIELD1",   i_colnames[1, RES_P_SFIELD1],  SZ_COLNAME)
	call strcpy( "SFIELD2",   i_colnames[1, RES_P_SFIELD1+1],  SZ_COLNAME)
	call strcpy( "SFIELD3",   i_colnames[1, RES_P_SFIELD1+2],  SZ_COLNAME)
	call strcpy( "SFIELD4",   i_colnames[1, RES_P_SFIELD1+3],  SZ_COLNAME)


	iferr {
	    nrows = tbhgti( tp, "gridrows")
	    ncols = tbhgti( tp, "gridcols")
	    call tbcfnd( tp, i_colnames, colptr, NUM_R_INFOCOLS)
	} then {
	    call tbtnam (tp, mess, SZ_RES_MESS)
	    call strcat( " is not a reseau file ", mess, SZ_RES_MESS)
	    call error (RES_ER_NOTRESFILE, mess)
	}
	do k = 1, NUM_R_INFOCOLS {
	    if (colptr[k] == NULL)
		call error (1, "column not found in reseau table")
	}

	# Create the structure for the reseau file
	call salloc( rp, LEN_RES_INFO + ncols*nrows*RES_COLSPERRES, 
						TY_LONG)
	# puts NULLs in
	call amovkl( NULL, Meml( rp ), ncols*nrows*RES_COLSPERRES )

	RES_FILE(rp)   = tp
	RES_NROWS(rp)  = nrows
	RES_NCOLS(rp)  = ncols
	RES_EXTENT(rp) = RES_COLSPERRES*nrows*ncols

	for (i = 1 ; i<= NUM_R_INFOCOLS ; i = i + 1){
	    RES_COLPTR(rp, i) = colptr[i]
	}

	tabcols = tbpsta( tp, TBL_NCOLS)
	for ( i=1; i<=tabcols; i=i+1) {
	    colpt = tbcnum( tp, i)
	    call tbcinf ( colpt, colnum, colname, colunits, 
				colfmt, coldatatype, collendata, collenfmt)
	    ip = 5
	    if ( strncmp( colname, "COL", 3) != 0 )
		next
	    else if ( ctoi( colname, ip, ival) == 0 )
		next
	    else if ( ival <1 || ival > ncols*nrows)
		next

	    if ( colname[4] == 'X' )
		RES_COLPTR( rp, NUM_R_INFOCOLS + ival ) = colpt
	    else if ( colname[4] == 'Y' )
		RES_COLPTR( rp, NUM_R_INFOCOLS +nrows*ncols+ival ) = colpt
	}

	# check we have enough columns defined
	for ( i=1; i<=RES_EXTENT(rp)+NUM_R_INFOCOLS; i=i+1) 
	    if ( RES_COLPTR( rp, i) == NULL )
		call error( 0, "reseau file is missing essential columns")


end		

# RS_CLOSE -- close reseau file

procedure rs_close ( rp ) 

pointer	rp	# i: pointer to res structure

errchk	tbtclo
begin

	if ( rp != NULL) {
	    if ( RES_FILE( rp ) != NULL ) {
	    	call tbtclo ( RES_FILE(rp) )
	    	RES_FILE( rp ) = NULL
	    }
	}
	
end

# RS_NEW -- create new reseau file. RS_OPEN is assumed to have been called
# The minimum number of entries is specified, together with the grid size.
# The required columns are created and the reseau file structure is
# initialised. 

procedure rs_new( num_entries, ncols, rows, cols, datatype, rp)

int	num_entries, ncols	# i: no. of entries and columns
int	rows, cols		# i: grid size of reseaux
int	datatype		# i: datatype of coord. info
pointer	rp			# i/o: reseau file structure pointer

#--

int	i
int	outcol, outrow, nampos
pointer	rc, rf, ru, rt, rl
pointer	tp, rpnew, sp

errchk	tbpset, tbtcre, tbcdef

begin
	tp = RES_FILE(rp)

	# decide how many row and columns to make
	outcol = NUM_R_INFOCOLS + RES_COLSPERRES*max(ncols, rows*cols) + 
		 NUM_R_EXTRACOLS
	outrow = RES_ROWSPERENTRY*(max(1, num_entries ) + NUM_R_EXTRAROWS)
	call tbpset( tp, TBL_ROWLEN, outcol)
	call tbpset( tp, TBL_ALLROWS, outrow)
	call tbpset( tp, TBL_MAXPAR, RES_USERPARS)
	RES_EXTENT(rp) = RES_COLSPERRES*rows*cols

	# create the table

	call tbtcre(tp)

	# expand reseau structure:	
	call salloc( rpnew, LEN_RES_INFO+RES_COLSPERRES*cols*rows, TY_LONG)
	# move the reseau info to the new structure
	call amovl(Meml[rp], Meml[rpnew], LEN_RES_INFO+RES_COLSPERRES*cols*rows)
	rp = rpnew

	# prepare to define all the reseau columns in one go
	call smark ( sp )
	call salloc( rc, outcol*(SZ_COLNAME+1),  TY_CHAR)
	call salloc( rf, outcol*(SZ_COLFMT+1),   TY_CHAR)
	call salloc( ru, outcol*(SZ_COLUNITS+1), TY_CHAR)
	call salloc( rt, outcol, TY_INT)
	call salloc( rl, outcol, TY_INT)


	call strcpy( "NAME", Memc[rc + (RES_P_ENTRY-1)*(SZ_COLNAME+1)], 
							SZ_COLNAME)
        Memi[rl + RES_P_ENTRY - 1] = 1
	Memi[rt + RES_P_ENTRY - 1] = -1*SZ_R_ENTRY

	call strcpy( "TRACKING",  
		Memc[rc + (RES_P_TRACKING-1)*(SZ_COLNAME+1)], SZ_COLNAME)
        Memi[rl + RES_P_TRACKING - 1] = 1
	Memi[rt + RES_P_TRACKING - 1] = -1*SZ_R_TRACKING

	call strcpy( "DATE", Memc[rc + (RES_P_DATE-1)*(SZ_COLNAME+1)],  
							SZ_COLNAME)
        Memi[rl + RES_P_DATE - 1] = 0
	Memi[rt + RES_P_DATE - 1] = TY_INT

	call strcpy( "NROWS", Memc[rc + (RES_P_NROWS-1)*(SZ_COLNAME+1)],  
							SZ_COLNAME)
        Memi[rl + RES_P_NROWS - 1] = 0
	Memi[rt + RES_P_NROWS - 1] = TY_INT

	call strcpy( "NCOLS", Memc[rc + (RES_P_NCOLS-1)*(SZ_COLNAME+1)],  
							SZ_COLNAME)
        Memi[rl + RES_P_NCOLS - 1] = 0
	Memi[rt + RES_P_NCOLS - 1] = TY_INT

	call strcpy( "SFIELD1",  Memc[rc + (RES_P_SFIELD1-1)*(SZ_COLNAME+1)],
							SZ_COLNAME)
        Memi[rl + RES_P_SFIELD1 - 1] = 0
	Memi[rt + RES_P_SFIELD1 - 1] = TY_REAL

	call strcpy( "SFIELD2",  Memc[rc + (RES_P_SFIELD2-1)*(SZ_COLNAME+1)],  
							SZ_COLNAME)
        Memi[rl + RES_P_SFIELD2 - 1] = 0
	Memi[rt + RES_P_SFIELD2 - 1] = TY_REAL

	call strcpy( "SFIELD3",  Memc[rc + (RES_P_SFIELD3-1)*(SZ_COLNAME+1)],
							SZ_COLNAME)
        Memi[rl + RES_P_SFIELD3 - 1] = 0
	Memi[rt + RES_P_SFIELD3 - 1] = TY_REAL

	call strcpy( "SFIELD4",  Memc[rc + (RES_P_SFIELD4-1)*(SZ_COLNAME+1)],  
							SZ_COLNAME)
        Memi[rl + RES_P_SFIELD4 - 1] = 0
	Memi[rt + RES_P_SFIELD4 - 1] = TY_REAL

	# Create the positional columns 
	for (i = 1; i<=cols*rows; i = i+1){

	    nampos = i + NUM_R_INFOCOLS - 1
	    # x coordinate column
	    call sprintf( Memc[rc+nampos*(SZ_COLNAME+1)], SZ_COLNAME, "COLX%u")
	    call pargi(i)
	    Memi[rl + nampos] = 0
	    Memi[rt + nampos] = RES_TYPE_TYPE

	    nampos = nampos + rows*cols
	    # y coordinate column
	    call sprintf( Memc[rc+nampos*(SZ_COLNAME+1)], SZ_COLNAME, "COLY%u")
	    call pargi(i)
	    Memi[rl + nampos] = 0
	    Memi[rt + nampos] = RES_TYPE_TYPE
	}

	# default units and format

	for (i = 1; i<=RES_COLSPERRES*cols*rows+NUM_R_INFOCOLS; i = i+1){
	    Memc[rf + (i-1)*(SZ_COLFMT+1)]   = EOS
	    Memc[ru + (i-1)*(SZ_COLUNITS+1)] = EOS
	}

	call tbcdef( tp, RES_COLPTR(rp, 1), 
			Memc[rc], Memc[ru], Memc[rf], 
			Memi[rt], Memi[rl], 
			RES_COLSPERRES*cols*rows+NUM_R_INFOCOLS)

	RES_FILE(rp) = tp
	RES_NROWS(rp) = rows
	RES_NCOLS(rp) = cols

	# Insert the SIZE info
	call tbhadi( tp, "gridrows", rows)
	call tbhadi( tp, "gridcols", cols)

	call sfree ( sp )

end                 

# RS_READR -- read reseau entry from reseau file.
# assumes rs_open has been called to get column pointers

int procedure rs_readr ( rp, rentp, ent_name, ent_num)

pointer	rp				# i: pointer to table
pointer	rentp				# i: pointer to res info structure
char	ent_name[SZ_R_PATTERN]		# i: entry name pattern requested
int	ent_num				# i: row requested - ignored if zero
#--
pointer	nullres				# temporary null flag pointer
pointer nullx				# temp null flag pointer for X
char	act_name[SZ_R_ENTRY]
int	nrows				# number of res grid rows
int	ncols				# number of res grid cols
char	tracking[SZ_R_TRACKING]		# associated tracking info 

int	i
int	numrows
int	next_num
int	act_num
int	icount
int	rs_num_rows()
bool	match
bool	nullflag[NUM_R_INFOCOLS]
bool	streq()
char	cinfo[SZ_R_CMAXSIZE, NUM_R_INFOCOLS]
int	num_to_read

errchk	tbrgtr, tbrgti, tbrgtt, rs_num_rows

begin

	act_num  = RES_F_NONEXTENTRY
	next_num = RES_F_NONEXTENTRY

	# just return if the entry row is negative
	if (ent_num < 0)		# "<=" changed to "<"; PEH 12/26/89
	    return (RES_F_NONEXTENTRY)

	numrows  = rs_num_rows(rp)

	# Check if row requested is beyond end of reseau entries
	if (ent_num > numrows || ent_num < 0)
	    call error(RES_ER_OVEND, "requested row beyond end of reseau file")

	# May be either access starting at row - id ent_num >0, 
	# or else by ent_name starting at beginning
	# PEH 1/3/90:  reorganize; don't take last entry if no match.

	if (ent_num <= 0) {		# find entry name

	    # If whitespace then take first entry.
	    if (ent_name[1] == EOS) {
		next_num = RES_F_FIRST
		act_num  = 1
		match    = true
	    } else {
		match    = false
		do icount = 1, numrows {
		    if (streq (ent_name,
			Memc[RES_ALL_NAMES(rp)+(icount-1)*(SZ_R_ENTRY+1)])) {

			match    = true
			next_num = RES_F_USEDWANTED
			act_num  = icount
			break
		    }
		}
	    }


	    if (!match) {
		call eprintf ("entry `%s' not found\n")
		    call pargstr (ent_name)
		call error (1, "")
	    }

	# If a specific row was requested then take that one 
	} else if (ent_num > 0) {
	    act_num  = ent_num
	    next_num = RES_F_USEDWANTED
	}

	# Now read the info about the entry wanted

	if (act_num != RES_F_NONEXTENTRY){
	    call tbrgtt( RES_FILE(rp), RES_COLPTR(rp, RES_P_ENTRY), 
			cinfo[1, RES_P_ENTRY], nullflag, SZ_R_ENTRY, 1, 
			(act_num-1)*RES_ROWSPERENTRY + 1)
	    call tbrgtt( RES_FILE(rp), RES_COLPTR(rp, RES_P_TRACKING), 
			cinfo[1, RES_P_TRACKING], nullflag, 
			SZ_R_TRACKING, 1, 
			(act_num-1)*RES_ROWSPERENTRY + 1)
	    call tbrgti( RES_FILE(rp), RES_COLPTR(rp, RES_P_DATE), 
			RES_ENT_DATE(rentp), nullflag, 1,
			(act_num-1)*RES_ROWSPERENTRY + 1)
	    call tbrgti( RES_FILE(rp), RES_COLPTR(rp, RES_P_NROWS), 
			RES_ENT_NROWS(rentp), nullflag, 1, 
			(act_num-1)*RES_ROWSPERENTRY + 1)
	    call tbrgti( RES_FILE(rp), RES_COLPTR(rp, RES_P_NCOLS), 
			RES_ENT_NCOLS(rentp), nullflag, 1, 
			(act_num-1)*RES_ROWSPERENTRY + 1)
	    call tbrgtr( RES_FILE(rp), 
			RES_COLPTR(rp, RES_P_SFIELD1), 
			RES_ENT_SFIELD1(rentp), nullflag, 
			NUM_R_RINFOCOLS, 
			(act_num-1)*RES_ROWSPERENTRY + 1)

	    call strcpy( cinfo[1, RES_P_ENTRY], act_name, SZ_R_ENTRY)
	    call strcpy( cinfo[1, RES_P_TRACKING], tracking, SZ_R_TRACKING)

	    nrows = RES_ENT_NROWS(rentp) 
	    ncols = RES_ENT_NCOLS(rentp)

	    # Now read position data.  Allocate memory for null flags for X.
	    num_to_read = nrows*ncols
	    call malloc (nullx, num_to_read, TY_BOOL)
	    nullres = RES_FPT(rentp)
	    call tbrgtr( RES_FILE(rp), RES_COLPTR(rp, NUM_R_INFOCOLS+1), 
			Memr[RES_XPT(rentp)], Memb[nullx], num_to_read, 
				(act_num-1)*RES_ROWSPERENTRY + 1)
	    call tbrgtr( RES_FILE(rp), 
			RES_COLPTR(rp, NUM_R_INFOCOLS+num_to_read+1), 
			Memr[RES_YPT(rentp)], Memb[nullres], num_to_read, 
				(act_num-1)*RES_ROWSPERENTRY + 1)

	    # Check for indef X or Y instead of just Y.  Added 3/20/92 by PEH.
	    do i = 0, num_to_read-1 {		# zero indexed
		if (Memb[nullx+i] || Memb[nullres+i]) {
		    Memr[RES_XPT(rentp)+i] = INDEF
		    Memr[RES_YPT(rentp)+i] = INDEF
		    Memb[RES_FPT(rentp)+i] = true
		}
	    }
	}

	if ( next_num != RES_F_NONEXTENTRY && act_num != RES_F_NONEXTENTRY )
	    next_num = act_num 
	else
	    next_num = RES_F_NONEXTENTRY

	call strcpy( act_name, RES_ENT_ENTRY(rentp),    SZ_R_ENTRY)
	call strcpy( tracking, RES_ENT_TRACKING(rentp), SZ_R_TRACKING)

	return ( next_num )

end


# WAS_PATTERN -- was the name given a template?

bool	procedure    was_pattern (name)

char	name[ARB]		# i: name to examine
#--
char	test_ch
int	stridx()

begin
	test_ch = CH_BOL
	if (stridx (test_ch, name) > 0)
	    return (true)

	test_ch = CH_NOT
	if (stridx (test_ch, name) > 0)
	    return (true)

	test_ch = CH_EOL
	if (stridx (test_ch, name) > 0)
	    return (true)

	test_ch = CH_ANY
	if (stridx (test_ch, name) > 0)
	    return (true)

	test_ch = CH_CLOSURE
	if (stridx (test_ch, name) > 0)
	    return (true)

	test_ch = CH_CCL
	if (stridx (test_ch, name) > 0)
	    return (true)

	test_ch = CH_CCLEND
	if (stridx (test_ch, name) > 0)
	    return (true)

	test_ch = CH_RANGE
	if (stridx (test_ch, name) > 0)
	    return (true)

	test_ch = CH_ESCAPE
	if (stridx (test_ch, name) > 0)
	    return (true)

	test_ch = CH_WHITESPACE
	if (stridx (test_ch, name) > 0)
	    return (true)

	test_ch = CH_IGNORECASE
	if (stridx (test_ch, name) > 0)
	    return (true)

	test_ch = CH_MATCHCASE
	if (stridx (test_ch, name) > 0)
	    return (true)

	return (false)
end

# RS_WRITE -- write reseau entry to reseau file.
# assumes rs_open has been called to get column pointers

procedure rs_write ( rp, rentp, row_num)

pointer	rp				# i: pointer to table
pointer	rentp				# i: pointer to res info structure
int	row_num				# i: entry to write, or next entry if 0
#--

int	tbl_row
int	numrow, num_to_write
int	rs_num_rows()
int	nrows, ncols

errchk	tbrptr, tbrpti, tbrptt, rs_num_rows

begin

	if (row_num == 0)
	    numrow = rs_num_rows(rp) + 1
	else
	    numrow = row_num

	nrows = RES_ENT_NROWS(rentp)
	ncols = RES_ENT_NCOLS(rentp)

	tbl_row = (numrow-1)*RES_ROWSPERENTRY + 1

#	Now write the info about the entry

	call tbrptt( RES_FILE(rp), RES_COLPTR(rp, RES_P_ENTRY), 
			RES_ENT_ENTRY(rentp), SZ_R_ENTRY, 1, tbl_row)

	call tbrptt( RES_FILE(rp), RES_COLPTR(rp, RES_P_TRACKING), 
			RES_ENT_TRACKING(rentp), 
			SZ_R_TRACKING, 1, tbl_row)

	call tbrpti( RES_FILE(rp), RES_COLPTR(rp, RES_P_DATE), 
			RES_ENT_DATE(rentp), 1, tbl_row)

	call tbrpti( RES_FILE(rp), RES_COLPTR(rp, RES_P_NROWS), 
			RES_ENT_NROWS(rentp), 1, tbl_row)
	call tbrpti( RES_FILE(rp), RES_COLPTR(rp, RES_P_NCOLS), 
			RES_ENT_NCOLS(rentp), 1, tbl_row)
	call tbrptr( RES_FILE(rp), 
			RES_COLPTR(rp, RES_P_SFIELD1), 
			RES_ENT_SFIELD1(rentp), 
			NUM_R_RINFOCOLS, tbl_row)

#	Now write position data

	num_to_write = RES_ENT_NROWS(rentp)*RES_ENT_NCOLS(rentp)

	call tbrptr( RES_FILE(rp), RES_COLPTR(rp, NUM_R_INFOCOLS+1), 
			Memr[RES_XPT(rentp)], num_to_write, tbl_row )
	call tbrptr( RES_FILE(rp), 
			RES_COLPTR(rp, NUM_R_INFOCOLS+1+num_to_write), 
			Memr[RES_YPT(rentp)], num_to_write, tbl_row )

end


# RS_INSERT -- insert entry into reseau file 
# If an input is given the entries are copied to the output before
# adding the new entry

procedure rs_insert(rpin, rp, respt)

pointer	rpin		# i: pointer to input reseau fie ( or NULL if none)
pointer	rp		# i: pointer to output res. file structure
pointer	respt		# i: reseau entry structure
#--

int	row_num
int	rs_readr(), junk, rs_num_rows()
pointer	rtemp, sp

errchk	rs_write, rs_readr, rs_num_rows, tbrgtt

begin

	# just write entry if there is no input file, 
	if (rpin != NULL && rpin != rp ) {
	    call smark(sp)
	    call rs_alloc( rp, RES_DATA_TYPE(respt), rtemp )
	    for (row_num = 1; row_num <= rs_num_rows(rpin);
					row_num = row_num + 1) {
		junk = rs_readr ( rpin, rtemp, EOS, row_num )
		call rs_write(rp, rtemp, row_num)
	    }
	    call sfree(sp)
	}

	call rs_write(rp, respt, 0)

end

# RS_CHECK -- check reseau size against given grid size

bool procedure rs_check( rp, nrows, ncols)

pointer	rp		# i: reseau file pointer
long	nrows		# i: number of rows in grid
long	ncols		# i: number of columns in grid

begin

	#see if total size matches nrows*ncols

	if (nrows*ncols == RES_NROWS(rp)*RES_NCOLS(rp) )
	    return(true)
	else
	    return(false)

end

# RS_NUM_ROWS -- return the number of reseau entries in the file

int procedure rs_num_rows(rp)

pointer	rp		# i: reseau file pointer
#--
errchk	tbpsta
int	tbpsta()

begin

	return(tbpsta(RES_FILE(rp), TBL_NROWS)/RES_ROWSPERENTRY)

end

# RS_T_UPDATE -- adds tracking info into reseau entry info block

procedure rs_t_update( rentp, tinfo, len)

pointer	rentp			# i: pointer to entry info block
char	tinfo[ARB]		# i: info to append
int	len			# i: max length of res info block
#--

long	clktime()

begin

	# add time stamp
	RES_ENT_DATE(rentp) = clktime(0)

	# initially just concatenate the info to the tracking info
	call strcat( tinfo, RES_ENT_TRACKING(rentp), len)

end

# XCOORD_TYPE -- given a file name, figure out if it is a text file or
# a table and then read the data into x and y arrays

long procedure xcoord_type( coordfile, entry, xcol, ycol, xcoord, ycoord )

char	coordfile[ARB]		# i: name of file with coordinates
char	entry[ARB]		# i: keyword for entry in reseau
char	xcol[ARB]		# i: keyword for xcol in table
char	ycol[ARB]		# i: keyword for ycol in table
pointer	xcoord			# o: pointer to x coord array
pointer	ycoord			# o: pointer to y coord array
#--

pointer	rs_open(), rs_template(), rlist, tbtopn(), rp, tp, respt
pointer	txcoord, tycoord, tcoord, xnull, ynull, nullpt, colptr[2]
int	rs_readr(), rs_list()
long	npts, nres, i
int	nrows, ncols, tbpsta()
char	ent_name[SZ_R_ENTRY], col[SZ_COLNAME, 2]
bool	streq()

errchk	rs_open, rs_template, rs_readr, rs_list
errchk	tbcgtr, xgfcoord

begin
	tp	= NULL
	rp	= NULL
	npts	= 0

	# If coordfile is the standard input, read it as a text file.
	if (streq (coordfile, "STDIN")) {

	    nrows = INDEFI
	    ncols = INDEFI
	    iferr( 
	      call xgfcoord( coordfile, INDEFI, INDEFI, nrows, ncols, tcoord) ) 
		call error (1, "error reading STDIN as coord file")
	    	
	    # transfer to 2 arrays
	    nres = nrows*ncols
	    npts = 0
	    call salloc ( xcoord, nrows*ncols, TY_LONG)
	    call salloc ( ycoord, nrows*ncols, TY_LONG)
	    for (i = 0 ; i <= nres - 1 ; i=i+1){
		if ( !IS_INDEF( Memr[tcoord+i] ) && 
		     !IS_INDEF( Memr[tcoord+i] )    ) {
		    Meml[xcoord + npts ] = nint ( Memr[tcoord+2*i] )
		    Meml[ycoord + npts ] = nint ( Memr[tcoord+2*i+1] )
		    npts = npts + 1
		}
	    }
	    return (npts)
	}

	# is it a RESEAU file?
	iferr ( rp = rs_open( coordfile, READ_ONLY, 0 ) )
	    rp = NULL

	if ( rp != NULL ) {
	    call clgstr( entry, ent_name, SZ_R_PATTERN)
	    call rs_alloc(  rp, TY_REAL, respt )
	    rlist = rs_template( rp, ent_name)
	    if (rs_readr( rp, respt, EOS, rs_list(rlist)) == RES_F_NONEXTENTRY )
	        call error ( 0, "could not read entry in reseau file")

	    # transfer non-null data to x,y long arrays
	    nres = RES_ENT_NROWS(respt)*RES_ENT_NCOLS(respt)
	    npts = 0
	    call salloc ( xcoord, nres, TY_LONG)
	    call salloc ( ycoord, nres, TY_LONG)
	    nullpt = RES_FPT(respt)
	    for (i = 0 ; i <= nres - 1 ; i=i+1){
		if ( !Memb[nullpt+i] ) {
		    Meml[xcoord + npts ] = nint ( Memr[RES_XPT(respt)+i] )
		    Meml[ycoord + npts ] = nint ( Memr[RES_YPT(respt)+i] )
		    npts = npts + 1
		}
	    }
	    call rs_close( rp )
	    return (npts)		
	} 

	# is this a TABLE 
	iferr ( tp = tbtopn( coordfile, READ_ONLY, 0) )
	    tp = NULL

	if ( tp != NULL) {
	    # get the column names to read from ordinary table
	    call clgstr( xcol, col[1, 1], SZ_COLNAME)
	    call clgstr( ycol, col[1, 2], SZ_COLNAME)
	    nrows = tbpsta( tp, TBL_NROWS)
	    call salloc ( txcoord, nrows, TY_REAL)
	    call salloc ( tycoord, nrows, TY_REAL)
	    call salloc ( xnull , nrows, TY_BOOL)
	    call salloc ( ynull , nrows, TY_BOOL)
	    call tbcfnd( tp, col, colptr, 2)
	    if (colptr[1] == NULL || colptr[2] == NULL)
		call error (1, "coord file is table, but column not found")
	    call tbcgtr( tp, colptr[1], Memr[txcoord], Memb[xnull], 1, nrows)
	    call tbcgtr( tp, colptr[2], Memr[tycoord], Memb[ynull], 1, nrows)
	    npts = 0
	    call salloc ( xcoord, nrows, TY_LONG)
	    call salloc ( ycoord, nrows, TY_LONG)
	    for ( i=0; i<=nrows-1; i=i+1) {
		if ( !Memb[xnull+i] && ! Memb[ynull+i] ) {
		    Meml[xcoord+npts] = nint (Memr[txcoord+i])
		    Meml[ycoord+npts] = nint (Memr[tycoord+i])
		    npts = npts + 1
		}
	    }
	    call tbtclo( tp )
	    return (npts)		
	}

	# Otherwise error.
	call error (1, "coord file is neither reseau, table nor text file")
end

# RS_ALLOC -- allocate space for reseau data

procedure rs_alloc( rp, datatype, respt)

pointer	rp		# i: res. structure pointer
int	datatype	# i: datatype for coordinate data
pointer	respt		# o: res. data pointer
#--

errchk	salloc

begin

	call salloc( respt, LEN_ENTRY_INFO, TY_STRUCT)
	call salloc( RES_XPT(respt), RES_EXTENT(rp), datatype)
	call salloc( RES_YPT(respt), RES_EXTENT(rp), datatype)
	call salloc( RES_FPT(respt), RES_EXTENT(rp), TY_BOOL)
	RES_DATA_TYPE(respt) = datatype

end
# RS_CNVTIME -- return reseau date as string
procedure rs_cnvtime( respt, datetime, max_size)

pointer	respt		# i: entry pointer
char	datetime[ARB]	# o: date/time string
int	max_size	# size of datetime
#--

long	timelong
begin

	timelong = RES_ENT_DATE( respt)
	call cnvtime( timelong, datetime, max_size)

end
# XGFCOORD -- read text file with x,y coordinate pairs - we may know how
# many coords to expect, or it can be left INDEF, in which case the
# number of rows is set to 1.

procedure xgfcoord( coordfile, rsize1, rsize2, nrows, ncols, xy )

char	coordfile[SZ_FNAME]	# i: coordfile name
int	rsize1, rsize2		# i: grid size if coord file given
pointer	nrows, ncols		# o: grid size
pointer	xy			# o: pointers to x,y real array
#--

pointer	open()
pointer	fdcoord			# text file pointer
pointer	xynew
bool	sized
char	line[SZ_LINE]
long	i_in_line, ip
real	dval
int	ctor()
bool	finished
long	init_size, saved_init_size
int	getline()

errchk	getline, ctor, open

begin


	if ( !IS_INDEFI (rsize1) || !IS_INDEFI (rsize2) ) {
	    if ( !IS_INDEFI (nrows) || !IS_INDEFI (ncols) ) {
		nrows = max( nrows, rsize2 )
		ncols = max( ncols, rsize2 )
	    } else {
		nrows = rsize2
		ncols = rsize1
	    }
	}

	if ( IS_INDEFI (nrows) || IS_INDEFI (ncols) ) {
	    sized = false
	    init_size = INIT_COORD_SIZE
	} else {
	    sized = true
	    init_size  = 2*nrows*ncols
	    saved_init_size = init_size
	}


	# read in the text data into buffer - guess initial size, and
	# expand buffer if required later
	call salloc( xy, init_size, TY_REAL)

	if (coordfile[1] != EOS ) {
	    i_in_line = 1
	    ip = 0
	    fdcoord = open(coordfile, READ_ONLY, TEXT_FILE)
	    finished = false
	    while ( getline(fdcoord, line) != EOF && !finished ) {
		i_in_line = 1
		while ( ctor(line, i_in_line, dval) >= 1 && line[1] != '#' ) {
		    if (ip > init_size) {
			call salloc( xynew, init_size + INC_COORD_SIZE, TY_REAL)
			call amovr( Memr[xy], Memr[xynew], init_size)
			init_size = init_size + INC_COORD_SIZE
			xy = xynew
		    }
		    Memr[xy + ip ] = dval
		    ip = ip + 1
		}
		if (sized && ip >= init_size ) 
		    finished = true
	    }
	    call close(fdcoord)

	    # set the size as 1*ip if indef sized
	    if (!sized || ip/2 != saved_init_size) {
		nrows= 1
		ncols = ip/2
	    }
	
	} 
end

# RS_TEMPLATE -- expand reseau entry templates
#
pointer procedure rs_template( rp, ent_name)

pointer	rp		# i: res struct pointer
char	ent_name[ARB]	# i: entry name template
#--

pointer	sel, msel, sp
int	icount, len, num_used
char	template[SZ_R_PATTERN]
int	rs_num_rows(), numrows

errchk	rs_num_rows, rs_texp

begin

	numrows = rs_num_rows(rp)
	len      = numrows + 1 
	num_used = 0
	call malloc( msel, len, TY_INT)
	
	call strcpy( ent_name, template, SZ_R_PATTERN)
	
	# temporary array for chosen points
	call smark ( sp )

	# match template
	call rs_texp( Memc[RES_ALL_NAMES(rp)], SZ_R_ENTRY, template,
			numrows, msel, len, num_used )

	call sfree( sp )

	# convert list 
	call salloc( sel, num_used + 1 + LEN_RES_TINDEX, TY_INT)
	RES_TEMPL_COUNT(sel) = LEN_RES_TINDEX
	for ( icount=1; icount <= num_used; icount=icount+1) {
	    Memi[sel + LEN_RES_TINDEX + icount -1] = Memi[msel + icount -1]
	}
	    
	call mfree ( msel, TY_INT )

	Memi[sel + num_used + LEN_RES_TINDEX] = RES_F_NONEXTENTRY
	
	return (sel)
end

# RS_LIST -- give next entry number to read

int procedure rs_list( rlist)
pointer	rlist		# i:   list pointer
#--

int	ip		# use count

begin

	ip = RES_TEMPL_COUNT(rlist) + 1
	if ( ip <= 0 )
	    call error(0, "illegal index in RES_LIST")

	# skip NULL entries
	while ( Memi[rlist +ip-1] == NULL )
	    ip = ip + 1

	RES_TEMPL_COUNT(rlist) = ip
	return ( Memi[rlist + ip - 1] )

end
# RS_TEXP -- Expand a entry template into an array of entry pointers
#
# Given an array of entry names and a entry name template,
# return an array of entry pointers. The size of the entry pointer array
# is given by nument and should be equal to the number of
# entrys in the table. 
#
# D. Giaretta	24-Jul-87	First Code

procedure rs_texp ( entnames, sz_ent, template, nument, msel, len, num_used)

char	entnames[sz_ent, ARB]	# i: array of entry names
int	sz_ent		# i: first dimension of entry name array
char	template[ARB]	# i: entry template
int	nument		# i: size of entry pointer array
pointer	msel		# i/o: pointer to array of entry pointers
int	len		# i/o: length of array of pointers
int	num_used	# i/o: number of entries selected
#--

bool	nometa		# true if pattern does not contain metacharacters
int	insert		# NO if template starts with negation character

int	fd_ptr		# pointer to stack of open list file descriptors
int	ic		# first non-white character in template

pointer fd_stack[MAX_STACK]
			# stack of file descriptors for open list files

pointer	sp, entpat, pattern, fd, tempent, tentpat

string	stkovflerr	"List files are nested too deeply, stack overflow"

int	strlen(), rstgetpat(), i, rststrmatch(), rstpatmatch(), numtmp
pointer	stropen(), open()

errchk	salloc, stropen, open, close
errchk	rstgetpat, rstmakpat, rststrmatch, rstpatmatch

begin

	call smark (sp)
	call salloc (entpat, SZ_FNAME, TY_CHAR)
	call salloc (pattern, SZ_FNAME, TY_CHAR)
	call salloc (tempent, nument, TY_INT)

	# Check the entry name template to find the first non-white character.

 	for (ic = 1; IS_WHITE (template[ic]); ic = ic + 1) 
		;
	
	if (template[ic] == EOS) {

####################need to decide about defaults here ##########
	    # If the template is blank, include all entrys in the array
	
	    num_used = nument
	    if ( num_used > len)
		call realloc( msel, num_used + len, TY_INT )
	    do i = 1, nument
		Memi[msel + i -1] = i

	    fd_ptr = 0

	} else {

	    # Open the entry name template as a file and push on
	    # the list file stack

	    fd_ptr = 1
	    fd_stack[1] =
		stropen (template[ic], strlen(template[ic]), READ_ONLY)

	}

	while (fd_ptr > 0) {

	    # Pop file descriptor off of the list file stack

	    fd = fd_stack[fd_ptr]
	    fd_ptr = fd_ptr - 1

	    # Loop over all entry patterns in the file

	    while (rstgetpat (fd, Memc[entpat], SZ_FNAME) > 0) {

		if (Memc[entpat] == '@') {

		    # If this pattern is a list file name, push the
		    # current descriptor on the stack and open the file

		    if (fd_ptr == MAX_STACK)
			call error (BOUNDS, stkovflerr)
		    fd_ptr = fd_ptr + 1
		    fd_stack[fd_ptr] = fd
		    fd = open (Memc[entpat+1], READ_ONLY, TEXT_FILE)

		} else {

		    # Otherwise, encode the pattern and search the table
		    # for matching entry names. To speed the search, use
		    # a special routine if the pattern does not include
		    # metacharacters

		    # If we are negating this group then we should
		    # remember this fact and jump to the next character
		    if ( Memc[entpat] == NEGCHAR ) {
			insert = NO
			tentpat = entpat + 1		    
		    } else {
			insert = YES
			tentpat = entpat
		    }

		    call rstmakpat (Memc[tentpat], Memc[pattern], SZ_FNAME,
				    nometa)
		    if (nometa) {
			numtmp = rststrmatch (entnames, sz_ent, Memc[pattern], 
					insert, nument, msel, len, num_used)
		    } else {
			numtmp = rstpatmatch (entnames, sz_ent, Memc[pattern], 
					insert, nument, msel, len, num_used)
		    }
		}
	    }
	    call close (fd)
	}

	call sfree (sp)
end

# RSTGETPAT -- Get next comma or whitespace delimeted pattern from file
#
# Copy characters into entpat until a field delimeter or the maximum number of
# characters is reached. The number of characters in entpat is returned as the
# value of the function, so the procedure which calls this one can test for
# the last field in the template.
#
# D. Giaretta	24-Jul-87	First Code

int procedure rstgetpat (fd, entpat, maxch) 

pointer	fd		# i: template file descriptor
char	entpat[ARB]	# o: pattern from entry name template
int	maxch		# i: maximum number of characters in field 
#--

char	ch		# next character from template
int	iq		# pointer to character in entpat

char	getc()

begin
	# Skip leading whitespace or commas

	ch = getc (fd, ch)
	while (IS_CNTRL(ch) || ch == BLANK || ch == DELIM)
	    ch = getc (fd, ch)

	# Copy characters to entpat. End when maxch is reached, or
	# when comma, whitespace, or EOF is found

	for (iq = 1; iq <= maxch; iq = iq + 1) {

	    if (IS_CNTRL(ch) || ch == BLANK || ch == DELIM || ch == EOF)
		break

	    entpat[iq] = ch
	    ch = getc (fd, ch)
	}
	entpat[iq] = EOS

	# If loop is terminated because of maxch, eat remaining characters
	# in field

	while (! IS_CNTRL(ch) && ch != BLANK && ch != DELIM && ch != EOF)
	    ch = getc (fd, ch)

	# Return number of characters in entpat

	return (iq-1)
end 

# RSTMAKPAT -- Encode the entry pattern
#
# Create the pattern used by the matching routines. Check for metacharacters
# (unescaped pattern matching characters) to see if the faster constant
# pattern routine can be used.
#
# D. Giaretta	24-Jul-87	First Code

procedure rstmakpat (entpat, pattern, maxch, nometa)

char	entpat[ARB]	# i: entry pattern string
char	pattern[ARB]	# o: Encoded pattern string
int	maxch		# i: Maximum length of encoded pattern string
bool	nometa		# o: True if no metacharacters in string
#--
int	ic, ip
pointer	sp, buffer, errtxt, ib

int	stridx(), strlen(), patmake()

string	patovflerr "entry pattern too long (%s)"
string	badpaterr  "entry pattern has bad syntax (%s)"

begin
	call smark (sp)
	call salloc (buffer, maxch, TY_CHAR)
	call salloc (errtxt, SZ_LINE, TY_CHAR)

	nometa = true
	ib = buffer

	# Copy the entry pattern to a temporary buffer

	for (ic = 1; entpat[ic] != EOS ; ic = ic + 1) {

	    # Copy escape sequences, but do not count as metacharacters

	    if (entpat[ic] == ESCAPE && entpat[ic+1] != EOS) {
		Memc[ib] = ESCAPE
		ib = ib + 1
		ic = ic + 1

	    # Covert '*' to '?*', count as a metacharacter

	    } else if (entpat[ic] == '*') {
		nometa = false
		Memc[ib] = '?'
		ib = ib + 1

	    # Check for the other metacharacters

	    } else if (stridx (entpat[ic], "[?{") > 0)
		nometa = false

	    Memc[ib] = entpat[ic]
	    ib = ib + 1
	}
	Memc[ib] = EOS

	# Check the buffer length against maximum pattern length

	if (strlen (Memc[buffer]) > maxch) {
	    call sprintf (Memc[errtxt], SZ_LINE, patovflerr)
	    call pargstr (entpat)
	    call error (BOUNDS, Memc[errtxt])
	}

	# If no metacharacters, strip escape sequences

	if (nometa) {
	    ip = 1
	    for (ib = buffer; Memc[ib] != EOS; ib = ib + 1) {
		if (Memc[ib] == ESCAPE && Memc[ib+1] != EOS)
		    ib = ib + 1
		pattern[ip] = Memc[ib]
		ip = ip + 1
	    }  
	    pattern[ip] = EOS

	# Otherwise, encode with patmake

	} else if (patmake (Memc[buffer], pattern, SZ_LINE) == ERR) {
	    call sprintf (Memc[errtxt], SZ_LINE, badpaterr)
	    call pargstr (entpat)
	    call error (SYNTAX, Memc[errtxt])
	}

	call sfree (sp)
end

# RSTSTRMATCH -- Add a entry pointer for a entry name to the array
# Used to match entry names when the entry pattern contains no
# metacharacters.
#
# D. Giaretta	24-Jul-87	First Code

int procedure rststrmatch (entnames, sz_ent, pattern, insert, nument, 
				msel, len, num_used)

char	entnames[sz_ent, ARB]	# i: array of entry names
int	sz_ent		# i: first dimension of entry name array
char	pattern[ARB]	# i: entry pattern
int	insert		# i: insert mode?
int	nument		# i: size of entry pointer array
pointer	msel		# i/o: pointer to array of entry pointers
int	len		# i/o: length of structure
int	num_used	# i/o: number of entries selected
#--

int	iptr, i
bool	streq()


begin

	# See if the entry name matches
	# Note we jump out when we match rather than searching for some
	# other matching entry
	# Either add to list, or else if not in insert mode we
	# delete any earlier occurence of selected entry
	for (iptr = 1; iptr<=nument; iptr = iptr + 1) {
	    if ( streq( pattern, entnames[1, iptr] ) ) {
		if ( insert == YES ) {
		    num_used = num_used + 1
	    	    if ( num_used > len) {
			len = len + num_used
			call realloc( msel, num_used + len, TY_INT )
		    }
		    Memi[msel + num_used -1] = iptr
		} else {
		    for ( i=1; i <= num_used; i=i+1) {
			if ( Memi[msel+i-1] == iptr ) 
			    Memi[msel+i-1] = NULL
		    }
		}
		break
	    }
	}
	return ( num_used )
end

# RSTPATMATCH -- Find entry pointers for entrys matching a pattern
#
# This routine is called when the entry pattern contains metacharacters.
#
# D. Giaretta	27-Jul-87	First Code

int procedure rstpatmatch (entnames, sz_ent, pattern, insert, nument, 
				msel, len, num_used )

char	entnames[sz_ent, ARB]	# i: array of entry names
int	sz_ent		# i: first dimension of entry name array
char	pattern[ARB]	# i: entry pattern
int	insert		# i: insert mode?
int	nument		# i: size of entry pointer array
pointer	msel		# i/o: pointer to array of entry pointers
int	len		# i/o: length of array
int	num_used	# i/o: number used in array
#--

int	ient, first_ch, last_ch, i
int	gpatmatch(), strlen()

errchk	gpatmatch

begin
	# Compare the entry pattern to each entry name in the table
	do ient = 1, nument 
	    # Check the entry name for a match
	    if (gpatmatch (entnames[1, ient], pattern, first_ch, last_ch) > 0) 
		if ((first_ch == 1) && (last_ch == strlen(entnames[1, ient]))) {
		    if ( insert == YES ) {
		    	num_used = num_used + 1
	    	    	if ( num_used > len) {
			    len = len + num_used
			    call realloc( msel, num_used + len, TY_INT )
		    	}
		    	Memi[msel + num_used -1] = ient
		    } else {
		    	for ( i=1; i <= num_used; i=i+1) {
			    if ( Memi[msel+i-1] == ient ) 
			    	Memi[msel+i-1] = NULL
		    	}
		    }
		}
	return ( num_used )
end

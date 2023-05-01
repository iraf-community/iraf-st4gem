include	<pattern.h>
include	<ctype.h>
include	<time.h>
include	<tbset.h>

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
RES_UTILS -- 

long 	= xcoord_type( coordfile, entry, xcol, ycol, xcoord, ycoord )
		given a file name, figure out if it is a text file or
		a table and then read the data into x and y arrays.
		returns number of coordinate pairs. If table file will 
		request the names "xcol" and "ycol" to read the x, y coords.
		Only non-INDEF coordinates are returned

.endhelp ____________________________________________________________


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

pointer	rlist, tbtopn(), tp, respt
pointer	txcoord, tycoord, tcoord, xnull, ynull, nullpt, colptr[2]
long	npts, nres, i
int	nrows, ncols, tbpsta()
char	col[SZ_COLNAME, 2]
bool	streq()

errchk	tbcgtr, xgfcoord

begin
	tp	= NULL
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


# XGFCOORD -- read text file with x,y coordinate pairs - we may know how
# many coords to expect, or it can be left INDEF, in which case the
# number of rows is set to 1.

# definitions for reading text files

define  INIT_COORD_SIZE 2000
define  INC_COORD_SIZE  1000

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

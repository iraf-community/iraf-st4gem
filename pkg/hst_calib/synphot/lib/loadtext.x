include	<tbset.h>

# LOADTEXT -- Load a table column column containing text.

procedure loadtext( file, colname, coldata, ndat )

char	file[ARB]		# i: Input file name
char	colname[ARB]		# i: Name of column
char	coldata[SZ_LINE, ARB]	# o: column strings
int	ndat			# i: max number of rows to read in

#--

int	nchar, nrow
int	nowhite(), tbpsta()
pointer	tp, ptr, nulflg, errmsg
pointer	tbtopn()

string	nocol	"Column %s not found in file %s"

begin

	call strlwr(file)
	nchar = nowhite( file, file, SZ_FNAME)

	tp = tbtopn(file, READ_ONLY, NULL)
	call tbcfnd(tp, colname, ptr, 1)

	nrow = tbpsta(tp, TBL_NROWS)
	ndat = min(ndat, nrow)
	call malloc(nulflg, ndat, TY_BOOL)
	call malloc( errmsg, SZ_LINE, TY_CHAR)

	if( ptr != 0) {
	   call tbcgtt(tp, ptr, coldata, Memb[nulflg], SZ_LINE, 1, ndat)
	} else {
	   call sprintf( Memc[errmsg], SZ_LINE, nocol )
	      call pargstr( colname )
	      call pargstr( file )
	   call error( 1, Memc[errmsg] )
	}

	call mfree( nulflg, TY_BOOL )
	call mfree( errmsg, TY_CHAR )

end

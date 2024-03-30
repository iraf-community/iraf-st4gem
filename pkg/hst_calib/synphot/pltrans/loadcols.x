include	<tbset.h>
include	"pltrans.h"

# LOADCOLS -- Read an STSDAS photometry table, sort by starname, then
#	      obsmode, create color differences for same starname, if
#	      needed.

procedure loadcols(fname, nrow, ncol, data, diff, xmode, xform, ymode, yform )

char	fname[ARB]		# i: Name of file to read
int	nrow			# io: Max rows on input; actual rows on output
int	ncol			# io: Max col on input; actual col on output
real	data[MAXROW, MAXCOL]	# o: data array
char	diff[ARB]		# o: true if data is difference (residual)
char	xmode[ARB]		# o: mode of x data
char	xform[ARB]		# o: form of x data
char	ymode[ARB]		# o: mode of y data
char	yform[ARB]		# o: form of y data

# Dec 1990 Dave Bazell  added iferr around calls to get xmode, ymode, xform,
# 	and yform header parameters.  gets from cl if not in header.

int	nr, ic
int	tbpsta(), strsearch()
pointer	tp, sp, flux1, flux2, error1, error2, nulflg, f1ptr, f2ptr
pointer	e1ptr, e2ptr
pointer	tbtopn()

string	nocol	"Column %s not found, setting to INDEFR.\n"

begin
	call strlwr( fname )
	if (strsearch(fname, "none") > 0 ) {
	   nrow = 0
	   ncol = 0
	   return
	}

	call smark( sp )
	
	# Open table and get column pointers
	tp = tbtopn( fname, READ_ONLY, NULL)
	call tbcfnd( tp, "FLUX1", f1ptr, 1)
	call tbcfnd( tp, "STATERROR1", e1ptr, 1)
	call tbcfnd( tp, "FLUX2", f2ptr, 1)
	call tbcfnd( tp, "STATERROR2", e2ptr, 1)

	nr = tbpsta( tp, TBL_NROWS )
	nrow = min( nr, nrow )

	call salloc( nulflg, nrow, TY_BOOL )
	call salloc( flux1, nrow, TY_REAL )
	call salloc( error1, nrow, TY_REAL )
	call salloc( flux2, nrow, TY_REAL )
	call salloc( error2, nrow, TY_REAL )

	# Read header parameters.  If no header parameter get values from
	# cl parameters
	iferr( call tbhgtt( tp, "XMODE", xmode, SZ_LINE) )
	   call clgstr( "xmode", xmode, SZ_LINE )
	iferr( call tbhgtt( tp, "XFORM", xform, SZ_FNAME) )
	   call clgstr( "xform", xform, SZ_FNAME )
	iferr( call tbhgtt( tp, "YMODE", ymode, SZ_LINE) )
	   call clgstr( "ymode", ymode, SZ_LINE )
	iferr( call tbhgtt( tp, "YFORM", yform, SZ_FNAME) )
	   call clgstr( "yform", yform, SZ_FNAME )
	iferr( call tbhgtt( tp, "DIFF", diff, 4 ) )
	   call strcpy( "notd",diff, 4 )

	# Read x-axis flux
	if ( f1ptr > 0 )
	   call tbcgtr( tp, f1ptr, Memr[flux1], Memb[nulflg], 1, nrow )
	else {
	   call printf(nocol)
	      call pargstr( "FLUX1" )
	   do ic = 1, nrow
	      Memr[flux1+ic-1] = INDEFR
	}

	# Read x-axis error
	if ( e1ptr > 0 )
	   call tbcgtr( tp, e1ptr, Memr[error1], Memb[nulflg], 1, nrow )
	else {
	   call printf(nocol)
	      call pargstr( "STATERROR1" )
	   do ic = 1, nrow
	      Memr[error1+ic-1] = INDEFR
	}

	# y-axis flux
	if ( f2ptr > 0 )
	   call tbcgtr( tp, f2ptr, Memr[flux2], Memb[nulflg], 1, nrow )
	else {
	   call printf(nocol)
	      call pargstr( "FLUX2" )
	   do ic = 1, nrow
	      Memr[flux2+ic-1] = INDEFR
	}

	# y-axis error
	if ( e2ptr > 0 )
	   call tbcgtr( tp, e2ptr, Memr[error2], Memb[nulflg], 1, nrow )
	else {
	   call printf(nocol)
	      call pargstr( "STATERROR2" )
	   do ic = 1, nrow
	      Memr[error2+ic-1] = INDEFR
	}

	do ic = 1, nrow {
	   data[ic,1] = Memr[flux1+ic-1]
	   data[ic,2] = Memr[error1+ic-1]
	   data[ic,3] = Memr[flux2+ic-1]
	   data[ic,4] = Memr[error2+ic-1]
	}

	call sfree( sp )
	call tbtclo( tp )

end

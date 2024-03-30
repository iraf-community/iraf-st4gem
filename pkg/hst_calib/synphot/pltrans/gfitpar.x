include	<pkg/dttext.h>

# GFITPAR -- Get the fit parameters output by GFIT1D.  

procedure gfitpar( cv, xmin, xmax )

pointer cv		# o: Pointer to the curfit structure
real	xmin		# o: Min value of x in fit
real	xmax		# o: Max value of x in fit

int	ic, dt, nrecs, ninfo, ncoeff, curve_type, recnum, nscan
int	dtmap(), dtgeti(), dtlocate(), scan()

pointer	coeff

real	dtgetr()

char	recname[SZ_FNAME], arg[SZ_FNAME], fitfile[SZ_FNAME]

string	cstr	"coeff%d"

begin

	# Get fitting parameter file name

	call clgstr( "fitfile", fitfile, SZ_FNAME )

	# Open database and get the record with fit parameters

	dt = dtmap( fitfile, READ_ONLY )
	nrecs = DT_NRECS( dt )

	if ( nrecs > 1 ) {

	   call eprintf( "%s has %d records.\n" )
	      call pargstr( fitfile )
	      call pargi( nrecs )

	   do ic = 1, nrecs {
	      call eprintf("%s\n")
	         call pargstr( DT_NAME(dt, ic) )
	   }

	   call eprintf("Enter desired record name:  ")
	   call flush(STDERR)
	   nscan = scan()
	      call gargstr( recname ,SZ_FNAME )

	} else

	   call strcpy( DT_NAME(dt, 1), recname, SZ_FNAME )

	# Get fit type, number of coefficients, min and max x values

	recnum = dtlocate( dt, recname )
	curve_type = dtgetr( dt, recnum, "coeff1" )
	ncoeff = dtgetr( dt, recnum, "coeff2" )
	xmin = dtgetr( dt, recnum, "coeff3" )
	xmax = dtgetr( dt, recnum, "coeff4" )

	# Allocate memory for the coefficient array

	call malloc( coeff, ncoeff, TY_REAL )

	# Get coefficients

	ninfo = dtgeti( dt, recnum, "ninfo" )
	do ic = 5, ninfo {
	   call sprintf( arg, SZ_FNAME, cstr )
	      call pargi( ic )
	   Memr[coeff+ic-5] = dtgetr( dt, recnum, arg)
	}

	# Set up the curfit structure

	call cvset( cv, curve_type, xmin, xmax, Memr[coeff], ncoeff )

	call mfree( coeff, TY_REAL )
end

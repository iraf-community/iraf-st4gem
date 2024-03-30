define	THRESH	8
define	MAXTEST	500

# PLLIMIT -- Determine reasonable plot limits

procedure pllimit( nplot, xplot, size, xsmall, xlarge)

int	nplot		# i: Number of data points
real	xplot[ARB]	# i: Array of data
real	size		# i: Units between plot limits
real	xsmall		# o: Lower limit for plot
real	xlarge		# o: Upper limit for plot

real	unit, scale, xmu, xsg, test
pointer	plot
int	ntest, itest, i, i2mid, isize, nvals
double	calc1, calc2
string	nularr	"Empty array in PLLIMIT."

#  adapted from 'plotlim' aug 1984 by keith horne at ioa
#  Dave Bazell Nov 1989  SPP version
#  Nov 1989  DB  Fix error in scaling mags by removing rescaling of data

begin

	if( nplot <= 0)
	   call error(1, nularr)

	unit = abs( size )

	call malloc( plot, nplot, TY_REAL)

	do i = 1, nplot
	   Memr[plot+i-1] = xplot[i]

	# Compute mean and sigma of plot data (use only MAXTEST points)

	calc1 = 0.d0
	calc2 = 0.d0
	ntest = max(nplot,MAXTEST)
	scale = real(nplot-1)/(ntest-1)
	nvals = 0

	# Calculate mean and std dev, rejecting INDEFRs
	do itest=0,ntest-1 {
	   i = nint( itest*scale ) + 1
	   if ( Memr[plot+i-1] != INDEFR ) {
	      calc1 = calc1 + Memr[plot+i-1]
	      calc2 = calc2 + Memr[plot+i-1]*Memr[plot+i-1]
	      nvals = nvals + 1
	   }
	}

	if (nvals > 0)
	   xmu = calc1/nvals

	xsg = calc2 - calc1*xmu

	if( nvals > 1)
	   xsg = xsg / (nvals-1)
	xsg = sqrt( max( xsg,unit/2. ) )

	# Find extrema, but reject thresh*sigma deviants

	xlarge = xmu
	xsmall = xmu
	do itest=0, ntest-1 {
	   i = nint( itest*scale ) + 1 
	   test = 0.
	   if( xsg > 0. && Memr[plot+i-1] != INDEFR ) {
	      test = (Memr[plot+i-1] - xmu) / xsg
	      if(abs(test) < THRESH ) {
	         xlarge = max(xlarge, Memr[plot+i-1])
	         xsmall = min(xsmall, Memr[plot+i-1])
	      }
	   }
        }

	# Round range to nearest integer and mid to nearest half-integer
	# multiple of unit

	if( unit > 0. ) {
	   isize = 1.1 * (xlarge-xsmall)/unit + 2.
	   i2mid = (xlarge+xsmall)/unit + 0.5
	   xlarge = 0.5 * unit * (i2mid+isize)
	   xsmall = 0.5 * unit * (i2mid-isize)
	}

	call mfree( plot, TY_REAL )
end

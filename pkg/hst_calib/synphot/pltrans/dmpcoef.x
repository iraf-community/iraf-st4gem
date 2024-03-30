include	"pltrans.h"

# DMPCOEF -- Print out the coefficients from the polynomial fit
# in PTRANS

procedure dmpcoef( filename, nterm, coeff, chisqr )

char	filename[ARB]		# i: output file name
int	nterm			# i: number of terms in fit
real	coeff[MAXFIT,ARB]	# i: coefficients of fits
real	chisqr[MAXFIT]		# i: chisqr values for fits

int	fd, nmin, nmax, ix, ic, jc
int	open()

begin

	fd = open( filename, NEW_FILE, TEXT_FILE )

	# Print the header
	call fprintf( fd, "    CHISQR  " )
	do ic = 1, nterm {
	   call fprintf( fd, "  Coeff%d " )
	      call pargi(ic)
	}
	call fprintf( fd, "\n")      

	nmin = max( 1, nterm - MAXFIT + 1 )
	nmax = min( nterm, MAXFIT)

	do ic = nmin, nmax {
	   ix = ic - nmin + 1
	   call fprintf( fd, "%11.5g" )
	      call pargr( chisqr[ix] )

	   do jc = 1, ic {
	      call fprintf( fd, "%10.5g" )
	         call pargr( coeff[ix,jc] )

	   }
	   call fprintf( fd, "\n" )
	}

	call close( fd )

end

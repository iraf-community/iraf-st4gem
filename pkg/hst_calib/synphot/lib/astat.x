define	TOSEC	1000.
define	IVP	($1 + ($2-1)*nvert)

# ASTAT -- Print status report for an amoeba iteration

procedure astat( time0, iter, ndat, npar, pvert, chi2, ibest, chi2avg, model )

int	time0		# i: Initial reference time
int	iter		# i: Current iteration
int	ndat		# i: Number of data points to fit
int	npar		# i: Number of parameters in fit
real	pvert[ARB]	# i: Values of vertices
real	chi2[ARB]	# i: Chisqr values for each vertex
int	ibest		# i: Index of best vertex in chisqr array
real	chi2avg		# i: Average chisqr over verticies
char	model[ARB]	# i: Model string

int	ndiv, nvert, iv, ip, time
long	dtime
long	cputime()

begin

	nvert = npar + 1

	# CPU usage
	time = cputime( dtime )

	# Print status report
	call printf("\nIter = %8d, CPU = %12g\n")
	   call pargi( iter )
	   call pargr( (time - time0)/TOSEC )
	ndiv = max ( 1, ndat - npar )
	call printf("NDAT = %8d, <chi^2>/%d = %12g, MIN(chi^2)/%d = %12g\n")
	   call pargi( ndat )
	   call pargi( ndiv )
	   call pargr( chi2avg/ndiv )
	   call pargi( ndiv )
	   call pargr( chi2[ibest]/ndiv )
	call printf("Model:  %s\n")
	   call pargstr( model )
	do iv = 1, nvert {
	   call printf("Chi^2 = %12g")
	      call pargr( chi2[iv]/ndiv )
	      do ip = 1, npar {
	         call printf(" %12g")
	             call pargr( pvert[IVP(iv,ip)] )
	      }
	      call printf("\n")
	}
	call flush( STDOUT )
end

include	"../lib/amoebapar.h"
define	EPS	1.e-5
define	IVP	($1 + ($2-1)*nvert)
define	MAXVERT	(MAXPAR+1)

# ASPECFIT -- Use amoeba downhill simplex method to fit spectral model to data

procedure aspecfit( model, data, npar, par, dpar, digits )

char	model[ARB]	# i: String specifying model to be fit
char	data[ARB]	# i: String specifying data to be fit
int	npar		# i: Number of parameters to fit
real	par[ARB]	# i: Initial parameter values
			# o: Final parameter values
real	dpar[ARB]	# i: Parameter offsets
int	digits[ARB]	# o: Digits of precision for each parameter

#--
int	niter, nvert, ip, iv, ibest, ndat, ic, time0, iter
int	strlen(), clgeti()
long	dtime
long	cputime()
real	delta, chi2[MAXVERT], chi2avg, pvert[MAXVERT*MAXPAR], ftol
real	parerr[MAXPAR]
real	clgetr()
bool	verbose, restart, first
data	first /true/
bool	clgetb(), resetpar()

int	specfunk()
extern	specfunk

# Include the common block for communication with specfunk
include "../lib/amoebafit.h"

begin

	if ( strlen(model) <= 0 )
	   call error( 1, "No model for amoeba" )
	if ( npar <= 0 )
	   call error( 1, "No parameters for amoeba fit")

	# Load common block for communication with specfunk
	call strcpy( model, fitmodel, SZ_LINE)
	call strcpy( data, fitdata, SZ_LINE)

	# Get cl parameters needed locally
	restart = clgetb( "restart" )
	niter = clgeti( "niter" )
	ftol = clgetr( "ftol" )

	# Report initial starting parameters
	call insertpar(par, model )
	call printf( "Start:  %s\n")
	   call pargstr( model )

	# Create initial simplex with centroid at starting parameters 
	# if one of the following is true
	# 1. restart flag is set to true or
	# 2. this is the first time through or
	# 3. parameter values have been changed
	if ( restart || first || resetpar( npar, par ) ) {
	   first = false
	   nvert = npar + 1
	   do ip = 1, npar {
	      delta = max( EPS * abs( par[ip] ) , abs( dpar[ip] ) )
	      do iv = 1, nvert
	         pvert[ IVP(iv,ip) ] = par[ip] - delta/nvert
	      pvert[ IVP(ip,ip) ] = pvert[ IVP(ip,ip) ] + delta
	   }
	}

	# Evaluate chisqr at each vertex
	do iv = 1, nvert {
	   do ip = 1, npar
	      par[ip] = pvert[ IVP(iv,ip) ]
	   call insertpar(par, model )
	   call specchi2( model, fitdata, chi2[iv], ndat )
	}


	# Initialize counters
	dtime = 0
	time0 = cputime( dtime )

	# Load the best fit parameter values and print status report
	call getbest( chi2, npar, pvert, par, chi2avg ,ibest, model )
	call astat( time0, 0, ndat, npar, pvert, chi2, ibest, chi2avg, model )

	verbose = clgetb( "verbose" )

	# Perform niter iterations.
	do ic = 1, niter {
	   iter = 1
	   call amoeba( pvert, chi2, nvert, npar, npar, ftol, specfunk, iter)
	   # Load the best fit parameters in to par[] array
	   call getbest( chi2, npar, pvert, par, chi2avg ,ibest, model )

	   # Report current statistics if verbose is set
	   if ( verbose )
	      call astat( time0, ic, ndat, npar, pvert, chi2, ibest, 
	                  chi2avg, model )
	}

	# Calculate error of each parameter from values of simplex vertices
	call simperr( par, npar, pvert, parerr )

	# Use errors to determine number of sig figs to use in output
	call etosf( npar, par, parerr, digits )

end

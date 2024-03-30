include	"amoebapar.h"
define	EPS	1.e-5
define	IVP	($1 + ($2-1)*nvert)
define	MAXVERT	(MAXPAR+1)

# AMOEBAFIT -- Use amoeba downhill simplex method to fit model to data

procedure amoebafit( model, data, npar, par, dpar )

char	model[ARB]	# i: String specifying model to be fit
char	data[ARB]	# i: String specifying data to be fit
int	npar		# i: Number of parameters to fit
real	par[ARB]	# i: Initial parameter values
			# o: Final parameter values
real	dpar[ARB]	# i: Parameter offsets

#--
int	niter, nvert, ip, iv, ibest, ndat, ic, time0
int	strlen(), clgeti()
long	dtime
long	cputime()
real	delta, chi2[MAXVERT], chi2avg, pvert[MAXVERT*MAXPAR], ftol
real	clgetr()
bool	verbose
bool	clgetb()

int	fitfunk()
extern	fitfunk

# Include the common block for communication with fitfunk
include "amoebafit.h"

begin

	if ( strlen(model) <= 0 )
	   call error( 1, "No model for amoeba" )
	if ( npar <= 0 )
	   call error( 1, "No parameters for amoeba fit")

	# Load common block for communication with fitfunk
	call strcpy( model, fitmodel, SZ_LINE)
	call strcpy( data, fitdata, SZ_LINE)

	# Get cl parameters needed locally
	niter = clgeti( "niter" )
	ftol = clgetr( "ftol" )

	# Report initial starting parameters
	call insertpar(par, model )
	call printf( "Start:  %s\n")
	   call pargstr( model )

	# Create initial simplex with centroid at starting parameters
	nvert = npar + 1
	do ip = 1, npar {
	   delta = max( EPS * abs( par[ip] ) , abs( dpar[ip] ) )
	   do iv = 1, nvert
	      pvert[ IVP(iv,ip) ] = par[ip] - delta/nvert
	   pvert[ IVP(ip,ip) ] = pvert[ IVP(ip,ip) ] + delta
	}

	# Evaluate chisqr at each vertex
	do iv = 1, nvert {
	   do ip = 1, npar
	      par[ip] = pvert[ IVP(iv,ip) ]
	   call insertpar(par, model )
	   call evalchi2( model, fitdata, chi2[iv], ndat )
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
	   call amoeba( pvert, chi2, nvert, npar, npar, ftol, fitfunk, 1)
	   # Load the best fit parameters in to par[] array
	   call getbest( chi2, npar, pvert, par, chi2avg ,ibest, model )

	   # Report current statistics if verbose is set
	   if ( verbose )
	      call astat( time0, ic, ndat, npar, pvert, chi2, ibest, 
	                  chi2avg, model )
	}
end

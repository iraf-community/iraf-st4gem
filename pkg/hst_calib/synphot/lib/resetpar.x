include	"amoebapar.h"
define	FRAC	1.e-6

# RESETPAR -- This function compares current parameter values with old
# values to see if any have changed.  If so returns true, otherwise it
# returns false.

bool procedure resetpar ( npar, par )

int	npar		# i: number of parameters
real	par[ARB]	# i: parameter value array

int	ic
real	oldpar[MAXPAR]
bool	reset

begin

	reset = false

	do ic = 1, npar {

	   # Check if parameter value has changed
	   if ( oldpar[ic] - par[ic] > FRAC * par[ic] )
	      reset = true

	   # Copy current parameter values to save array
	   oldpar[ic] = par[ic]

	}
	return reset
end


# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<error.h>
include <mach.h>
include	"../at.h"
include	"../atom.h"
include	"../neberr.h"

#--------------------------------------------------------------------10 Jul 97--
.help temden.x May96 nebular/temden
.ih
NAME
         temden - Employs geometric bisection to find T_e or N_e 
set_iter_limits - Set upper/lower bounds for T_e/N_e calculation 
prt_verb_header - Print header for verbose output.
    prt_verbose - Print verbose output
.endhelp
#-------------------------------------------------------------------------------
#  TEMDEN -	Employs an iterative technique (geometric bisection) to find 
#		T_e or N_e.  Used only when T_e or N_e is unknown.  

int procedure temden (at, den, tem, diag_type, robs, diag_expr, errx, verbos)

#  Calling arguments:
pointer	at			# atomic data structure
real	den			# electron density
real	tem			# electron temerature
int	diag_type		# solve for DENSITY or TEMPERATURE
real	robs			# measured (observed) ratio from spectrum
char	diag_expr[ARB]		# expression for diagnostic ratio
real	errx			# fractional tolerance in X
bool	verbos			# print verbose output?

#  Local variables:
int	iter			# iterative index used to find TEM or DEN
real	rcal			# characteristic ratio calculated in SOLVE
real	rdev1, rdev2		# fractional deviations of obs/calc ratios
real	rlo, rhi		# ratios calculated at XLO,XHI
real	x			# current best guess for TEM or DEN
real	xlo, xhi		# current range on iterated TEM or DEN

#  Function used:
real	j_ratio()		# ratio of calculated line emissivities

include	"../flerr.com"

define	N_LVLS		  5	# No. atomic energy levels used
define	MAX_ITER	100	# 
define	DEVIATION    ( abs (0.5 * ($1-$2) / ($1+$2)) )

errchk	j_ratio, set_iter_limits

begin
	# Set iteration limits.
	call set_iter_limits (at, diag_type, xlo, xhi)

	if (verbos) 
	    call prt_verb_header (diag_type)

	# Use geometric bisection to match calculated to observed ratio.
	if (IS_INDEFR(robs))
	    return (DIAG_INVALID)

	do iter = 1, MAX_ITER {
	    if (iter == 1)
		x = xlo
	    else if (iter == 2)
		x = xhi
	    else
		x = sqrt (xlo * xhi)

	    if (diag_type == TEMPERATURE) 
		tem = x

	    else if (diag_type == DENSITY) 
		den = x

	    iferr (call solve (at, den, tem))
		return (DIAG_FAILED)

	    rcal = j_ratio (at, diag_expr)
	    if (IS_INDEFR(rcal))
		return (DIAG_INVALID)

	    # Output intermediate iterations.
	    if (verbos) 
	    	call prt_verbose (diag_type, iter, den, tem, xlo, xhi, rcal, robs)

	    # Test for convergence.
	    if (iter == 1) {
		rlo = rcal
		next

	    } else if (iter == 2) {
		rhi = rcal

		# Return error if the input ratio is out of bounds.
	    	rdev1 = (rcal - robs) / (rlo - robs)
	    	if (rdev1 > 0.) {

		    # Write ratio & ranges to error common block.
		    val1 = robs
		    val2 = real (AT_ATOM(at))
		    val3 = min (rlo, rcal)
		    val4 = max (rlo, rcal)
		    return (LR_OUT_OF_BOUNDS)
		}
		next
	    }

	    # Changed convergence condition to account for zero deviation 
	    # of the observed ratio from that calculated. 
	    rdev1 = (rcal - robs) / (rlo - robs)
	    rdev2 = (rcal - robs) / (rhi - robs)
	    if (DEVIATION (xhi, xlo) < errx || abs (rdev2) <= EPSILONR) 

		# Converged: all done
		return (OK)

	    else if ( rdev2 < 0. ) {
		rlo = rcal
		xlo = x

	    } else if ( rdev2 > 0. ) {
		rhi = rcal
		xhi = x
	    }
	}

	# Did not converge after max_iterations. 
	return (NOT_CONVERGENT)

end


#-------------------------------------------------------------------------------
#  SET_ITER_LIMITS - Set upper/lower bounds for temperature/density calculation.

procedure set_iter_limits (at, diag_type, x_lo, x_hi)

#  Arguments:
pointer	at		# I: atomic data object
int	diag_type	# I: solve for DENSITY or TEMPERATURE
real	x_lo, x_hi	# O: lower, upper limits on TEM or DEN

#  Declarations:
int	ion		# spectrum code

begin
	# Temperature limits set by input collision strenth tables. 
	if (diag_type == TEMPERATURE) {
	    x_lo = AT_TE_MIN(at)
	    x_hi = AT_TE_MAX(at)

	# Density limits: 1. < dens < 1.e+8 /cm^3. 
	} else if (diag_type == DENSITY) {
	    x_lo = 1.
	    x_hi = 1.e8
	}

	ion = AT_ION(at)
 
	# Special cases:
	switch (AT_ATOM(at)) {
	case CARBON:
	    if (ion == 1) {
		if (diag_type == DENSITY) 
	    	    x_hi = 1.e5
	    }

	case NEON:
	    if (ion == 2 || ion == 4) {
		if (diag_type == TEMPERATURE) 
	    	    x_lo = 2000.
	    }

	case MAGNESIUM:
	    if (ion == 4 || ion == 6) {
		if (diag_type == TEMPERATURE) 
	    	    x_lo = 3000.
	    }

	case POTASSIUM:
	    if (ion == 3) {
		if (diag_type == TEMPERATURE) 
	    	    x_lo = 2000.
	    }

	case CALCIUM:
	    if (ion == 4 || ion == 6) {
		if (diag_type == TEMPERATURE) 
	    	    x_lo = 2000.
	    }
	}
end


#-------------------------------------------------------------------------------
#  PRT_VERB_HEADER - Print header for verbose output.

procedure prt_verb_header (diag_type)

#  Arguments:
int	diag_type	# solve for DENSITY or TEMPERATURE

begin
	call printf (
		"Iter#  %11s      XLO    %11s      XHI    Fract. Dev.\n")

	if (diag_type == TEMPERATURE) {
	    call pargstr ("Temperature")
	    call pargstr ("  Density  ")

	} else if (diag_type == DENSITY) {
	    call pargstr ("  Density  ")
	    call pargstr ("Temperature")
	}
end


#-------------------------------------------------------------------------------
#  PRT_VERBOSE - Print verbose output. 

procedure prt_verbose (diag_type, iter, den, tem, xlo, xhi, rcal, robs)

#  Arguments:
int	diag_type	# solve for DENSITY or TEMPERATURE
int	iter		# iteratation number
real	den		# electron density
real	tem		# electron temerature
real	xlo, xhi	# current range on iterated TEM or DEN
real	rcal		# characteristic ratio calculated in SOLVE
real	robs		# measured (observed) ratio from spectrum

begin
	call printf (" %3d   %10.5e  %10.5e  %10.5e  %10.5e   %9.4g \n")
	    call pargi (iter)
	if (diag_type == DENSITY) {
	    call pargr (den)
	    call pargr (xlo)
	    call pargr (tem)
	} else if (diag_type == TEMPERATURE) {
	    call pargr (tem)
	    call pargr (xlo)
	    call pargr (den)
    	}
	    call pargr (xhi)
	    call pargr (abs(rcal-robs)/(rcal+robs))
end



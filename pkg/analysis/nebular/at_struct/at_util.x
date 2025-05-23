# Copyright restrictions apply - see st4gem$copyright.stsdas

include <mach.h>
include "../at.h"

#--------------------------------------------------------------------13 Feb 98--
.help at_util.x Feb98 nebular/at_struct
.ih
NAME
.nf
      atg_wave - Compute the wavelength for an atomic transition
     atg_emiss - Compute the emissivity for an atomic transition
atg_valid_lvls - Are the upper & lower levels for a transition valid?
.fi
.ih
DESCRIPTION
These routines set or get various attributes for an atomic 
data structure. 
.endhelp
#-------------------------------------------------------------------------------
#  ATG_WAVE --	Compute the wavelength for an atomic transition. Returns 
#		INDEFD for invalid indicies. 

double procedure atg_wave (at, upper, lower)

#  Calling arguments:
pointer	at			# I: atomic data object
int	upper, lower		# I: energy level indexes
double	wave			# O: wavelength of transition in Angstroms

# Declaration:
bool	atg_valid_lvls()	# Validate upper/lower levels for a transition
int	n_lvls			# number of levels for atom

#  Memory management:
define	E_LVL		Memd[AT_E_TR($1)+($2-1)]

begin
	# Ensure that the upper/lower values are within the array bounds.
	n_lvls = AT_NLVL(at)
	if (!atg_valid_lvls (upper, lower, n_lvls))
	    return (INDEFD)

	wave = 1.d+0 / abs (E_LVL(at, upper) - E_LVL(at, lower))
	return (wave)
end


#-------------------------------------------------------------------------------
#  ATG_EMISS --	Compute the emissivity for an atomic transition. Returns 
#		INDEFD for invalid indicies. 

double procedure atg_emiss (at, upper, lower, density)

#  Calling arguments:
pointer	at			# I: atomic data object
int	upper, lower		# I: energy level indexes
real	density			# I: electron density
double	emissivity		# O: emissivity of transition 

# Declarations:
bool	atg_valid_lvls ()	# Validate upper/lower levels for a transition
int	n_lvls			# number of levels for atom

#  Memory management:
define	EM		Memd[AT_EMISS(at)+($1-1)+(($2-1)*AT_NLVL(at))]
define	MIN_DOUBLE	log10 (1. / MAX_DOUBLE)

begin
	# Ensure that the upper/lower values are within the array bounds.
	n_lvls = AT_NLVL(at)
	if (!atg_valid_lvls (upper, lower, n_lvls))
	    return (INDEFD)

	if (log10(EM(upper,lower)) - log10(double(density)) > MIN_DOUBLE)
	    emissivity = EM(upper,lower) / double (density)
	else
	    return (INDEFD)

	return (emissivity)
end


#-------------------------------------------------------------------------------
#  ATG_VALID_LVLS --	Validate the upper & lower levels for a transition. 

bool procedure atg_valid_lvls (upper, lower, n_lvls)

#  Calling arguments:
int	n_lvls			# I: number of levels for atom
int	upper, lower		# I: energy level indexes

begin
	# Ensure that the upper/lower values are within the array bounds.
	if (upper > n_lvls || upper < 2)
	    return (false)

	if (lower < 1 || lower > n_lvls - 1)
	    return (false)

	if (upper == lower)
	    return (false)

	return (true)
end



# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<error.h>
include	"../atom.h"
include	"../fivel.h"
include	"../zones.h"

define	MAX_AVG		5	# Max no. diagnostics to average
define	TOLER		1.e-5	# tolerance for computation
define	VERB		false	# verbose option for calc. iterations

#--------------------------------------------------------------------20 Feb 97--
.help zones.x May96 nebular/fivel
.ih
NAME
    z_init - Initialize ZONES data structure.  
     zone1 - Calculate T_e, N_e, from O, N, & S for low ionization zone 
     zone2 - Calculate T_e, N_e, from O, Ne, Ar, S, Cl for medium ionization zone 
     zone3 - Calculate T_e, N_e, from C, Ne, & Ar for high ionization zone 
zcalc_dens - Calculates pre-defined N_e from diagnostic 
zcalc_temp - Calculates pre-defined T_e from diagnostic 
.ih
DESCRIPTION
.endhelp
#-------------------------------------------------------------------------------
#  Z_INIT -	Initialize ZONES data structure.  

procedure z_init (zn)

#  Calling argument:
pointer	zn			# zones structure

begin
	Ne_Ciii(zn)  = INDEFR
	Ne_Oii(zn)   = INDEFR
	Ne_NEiv(zn)  = INDEFR
	Ne_Sii(zn)   = INDEFR
	Ne_CLiii(zn) = INDEFR
	Ne_ARiv(zn)  = INDEFR
	Te_Nii(zn)   = INDEFR
	Te_Oii(zn)   = INDEFR
	Te_Oiii(zn)  = INDEFR
	Te_NEiii(zn) = INDEFR
#	Te_NEiv(zn)  = INDEFR
	Te_NEv(zn)   = INDEFR
	Te_Sii(zn)   = INDEFR
	Te_Siii(zn)  = INDEFR
	Te_ARiii(zn) = INDEFR
	Te_ARiv(zn)  = INDEFR
	Te_ARv(zn)   = INDEFR

	Ne_Low(zn) = INDEFR
	Ne_Med(zn) = INDEFR
	Ne_Hi(zn)  = INDEFR
	Te_Low(zn) = INDEFR
	Te_Med(zn) = INDEFR
	Te_Hi(zn)  = INDEFR
end


#-------------------------------------------------------------------------------
#  ZONE1 -	Calculate T_e, N_e, from O, N, and S, for low ionization zone 
#		of nebula and update the output table.  

procedure zone1 (zn, ilr)

#  Calling arguments
pointer	zn			# zones structure
pointer	ilr			# line ratio structure

#  Local variables:
pointer	atl			# atomic data object list
real	density			# electron density
real	tempr			# electron temperature
real	wt[MAX_AVG]		# relative weights
real	val[MAX_AVG]		# working variables

#  Functions used:
real	zcalc_dens()		# Calculate N_e from diagnostic
real	zcalc_temp()		# Calculate T_e from diagnostic
real	wt_avg()		# Weighted avg of array, excluding INDEFs

errchk	zcalc_dens, zcalc_temp

begin
	atl = Z_ATL(zn)
	call amovkr (1., wt, MAX_AVG)

	# Calculate primary, low I.P. diagnostics.  
	# Begin by assuming N_e = 1000 and derive T_e[Nii].  

	density = 1000.

	# Calculate T_e from Nii and Oiii.
	val[1] = zcalc_temp (atl, NITROGEN, 1, density, RT_N2(ilr))
	val[2] = zcalc_temp (atl, OXYGEN,   2, density, RT_O3(ilr))

	# Average the two temperatures for first approximation to Te_low 
	wt[1] = 1.
	wt[2] = 2.
	tempr = wt_avg (val, wt, 2)
	if (IS_INDEF (tempr))
	    tempr = 10000.

	# Use the approximation to Te_low to calculate N_e. 
	# Calculate densities from Oii and Sii. 
	val[1] = zcalc_dens (atl, OXYGEN, 1, tempr, RN_O2(ilr))
	val[2] = zcalc_dens (atl, SULFUR, 1, tempr, RN_S2(ilr))

	# Use the average of Ne_[Oii] and Ne_[Sii] for Ne_Low, if available; 
	# otherwise, adopt Ne_Low = 1000 / cm^3.
	call amovkr (1., wt, MAX_AVG)
	density = wt_avg (val, wt, 2)
	if (IS_INDEFR (density)) 
	    density = 1000.

	# Now re-calculate T_Nii and T_Oiii using N_Low.  
	tempr       = zcalc_temp (atl, NITROGEN, 1, density, RT_N2(ilr))
	Te_Oiii(zn) = zcalc_temp (atl, OXYGEN,   2, density, RT_O3(ilr))

	# If tempr unavailable, use T_Oiii or 10,000 K.  
	if ( IS_INDEFR (tempr) ) {
	    if ( !IS_INDEFR (Te_Oiii(zn)) )
	    	tempr = Te_Oiii(zn)
	    else
	    	tempr = 10000.
	}

	# Finally, use tempr (which is T_Nii if available) to calculate 
	# Ne_Low = average of Ne_Oii and Ne_Sii
	val[1] = zcalc_dens (atl, OXYGEN, 1, tempr, RN_O2(ilr))
	val[2] = zcalc_dens (atl, SULFUR, 1, tempr, RN_S2(ilr))

	Ne_Oii(zn) = val[1]
	Ne_Sii(zn) = val[2]

	# Use the average of Ne_[Oii] and Ne_[Sii] for Ne_Low, if available 
	call amovkr (1., wt, MAX_AVG)
	Ne_Low(zn) = wt_avg (val, wt, 2)

	# Calculate Te_Oii and Te_Sii with improved density estimate. 
	if ( !IS_INDEF (Ne_Low(zn)) ) 
	    density = Ne_Low(zn)

	val[1] = zcalc_temp (atl, NITROGEN, 1, density, RT_N2(ilr))
	val[2] = zcalc_temp (atl, OXYGEN,   1, density, RT_O2(ilr))
	val[3] = zcalc_temp (atl, SULFUR,   1, density, RT_S2(ilr))

	Te_Nii(zn) = val[1]
	Te_Oii(zn) = val[2]
	Te_Sii(zn) = val[3]

	# Finally, average temperatures for T_low. 
	call amovkr (1., wt, MAX_AVG)
	wt[1] = 5
	Te_Low(zn) = wt_avg (val, wt, 3)

	end


#-------------------------------------------------------------------------------
#  ZONE2 -	Calculate T_e, N_e, from O, Ne, Ar, S, & Cl, for intermediate 
#		ionization zone of nebula.  

procedure zone2 (zn, ilr)

#  Calling arguments
pointer	zn			# zones structure
pointer	ilr			# line ratio structure

#  Local variables
pointer	atl			# atomic data object list
real	density			# electron density
real	tempr			# electron temperature
real	wt[MAX_AVG]		# relative weights
real	val[MAX_AVG]		# working variables

#  Functions used:
real	wt_avg()		# Weighted avg of array, excluding INDEFs
real	zcalc_dens()		# Calculate N_e from diagnostic
real	zcalc_temp()		# Calculate T_e from diagnostic

errchk	zcalc_dens, zcalc_temp

begin
	atl = Z_ATL(zn)
	call amovkr (1., wt, MAX_AVG)

	# Calculate primary, medium I.P. diagnostics. 
	# Begin by assuming N_e = 1000. 
	density = 1000.

	# Calculate estimate for T_e from T_e[O iii].  
	tempr = zcalc_temp (atl, OXYGEN, 2, density, RT_O3(ilr))

	if ( IS_INDEFR (tempr) ) 
	    tempr = 10000.

	# Use approximate T_Med to calculate Ne_[CLiii], Ne_[ARiv],  Ne_Ciii].  
	val[1] = zcalc_dens (atl, CHLORINE, 2, tempr, RN_CL3(ilr))
	val[2] = zcalc_dens (atl, ARGON,    3, tempr, RN_AR4(ilr))
	val[3] = zcalc_dens (atl, CARBON,   2, tempr, RN_C3(ilr))

	# Use the average of Ne_[Cl iii] and Ne_[Ar iv], if available, as N_Med. 
	# Note: relative weights revised 28-Apr-94
	call amovkr (1., wt, MAX_AVG)
	wt[1] = 2.0
	density = wt_avg (val, wt, 3)

	# Now re-calculate T_e from [O iii], [NE iii], [AR iii], & [S iii] 
	# using average density, if available.  
	if ( IS_INDEFR (density) ) 
	    density = 1000.

	val[1] = zcalc_temp (atl, OXYGEN, 2, density, RT_O3(ilr))
	val[2] = zcalc_temp (atl, NEON,   2, density, RT_NE3(ilr))
	val[3] = zcalc_temp (atl, ARGON,  2, density, RT_AR3(ilr))
	val[4] = zcalc_temp (atl, ARGON,  3, density, RT_AR4(ilr))
	val[5] = zcalc_temp (atl, SULFUR, 2, density, RT_S3(ilr))

	Te_Oiii(zn)  = val[1]
	Te_NEiii(zn) = val[2]
	Te_ARiii(zn) = val[3]
	Te_ARiv(zn)  = val[4]
	Te_Siii(zn)  = val[5]

	# Calculate average temperature for T_Med
	# Note: relative weights revised 12-Jan-96
	call amovkr (1., wt, MAX_AVG)
	wt[1] = 4.
	wt[2] = 2.
	wt[3] = 2.
	Te_Med(zn) = wt_avg (val, wt, 5)

	# Now improve the C iii], [CL iii] and [AR iv] density calculations.  
	if ( IS_INDEFR (Te_Med(zn)) ) 
	    tempr = 10000.
	else
	    tempr = Te_Med(zn)

	val[1] = zcalc_dens (atl, CHLORINE, 2, tempr, RN_CL3(ilr))
	val[2] = zcalc_dens (atl, ARGON,    3, tempr, RN_AR4(ilr))
	val[3] = zcalc_dens (atl, CARBON,   2, tempr, RN_C3(ilr))

	Ne_CLiii(zn) = val[1]
	Ne_ARiv(zn)  = val[2]
	Ne_Ciii(zn)  = val[3]

	# Average Ne_[Cl iii] and Ne_[Ar iv], if available, for Ne_Med. 
	call amovkr (1., wt, MAX_AVG)
	wt[3] = 0.5
	Ne_Med(zn) = wt_avg (val, wt, 3)

	end


#-------------------------------------------------------------------------------
#  ZONE3 -	Calculate T_e, N_e from C, Ne, & Ar for high ionization zone 
#		of nebula.  

procedure zone3 (zn, ilr)

#  Calling arguments
pointer	zn			# zones structure
pointer	ilr			# line ratio structure

#  Local variables
pointer	atl			# atomic data object list
real	density			# electron density
real	wt[MAX_AVG]		# relative weights
real	val[MAX_AVG]		# working variables

#  Functions used:
real	wt_avg()		# Weighted avg of array, excluding INDEFs
real	zcalc_dens()		# Calculate N_e from diagnostic
real	zcalc_temp()		# Calculate T_e from diagnostic

errchk	zcalc_dens, zcalc_temp

begin
	atl = Z_ATL(zn)
	call amovkr (1., wt, MAX_AVG)

	# Calculate high I.P. diagnostics.  
	# Assume T_e = Te_[O iii]. 
	if ( IS_INDEFR (Te_Oiii(zn)) )
	    Te_Hi(zn) = 10000.
	else
	    Te_Hi(zn) = Te_Oiii(zn)

	# Calculate Ne_[NEiv]. 
	Ne_NEiv(zn) = zcalc_dens (atl, NEON, 3, Te_Hi(zn), RN_NE4(ilr))

	Ne_Hi(zn) = Ne_NEiv(zn)
	if (IS_INDEFR (Ne_Hi(zn)) )
	    density = 1000.
	else
	    density = Ne_Hi(zn)

	# Use approximate Ne_Hi to calculate T_e from [Ar v], [Ne v], & [Ne iv]. 
	# Note: Te_[NEiv] not yet available.
	val[1] = zcalc_temp (atl, ARGON, 4, density, RT_AR5(ilr))
	val[2] = zcalc_temp (atl, NEON,  4, density, RT_NE5(ilr))
#	val[3] = zcalc_temp (atl, NEON,  3, density, RT_NE4(ilr))

	# Average the two temperatures for an approximation to T_Hi. 
	Te_ARv(zn)  = val[1]
	Te_NEv(zn)  = val[2]
#	Te_NEiv(zn) = val[3]
	call amovkr (1., wt, MAX_AVG)

	Te_Hi(zn) = wt_avg (val, wt, 2)

	# Now improve the [NEiv] density calculation. 
	Ne_NEiv(zn) = zcalc_dens (atl, NEON, 3, Te_Hi(zn), RN_NE4(ilr))

	Ne_Hi(zn) = Ne_NEiv(zn)
	call flush (STDOUT)
end


#-------------------------------------------------------------------------------
#  ZCALC_DENS -	Calculates pre-defined N_e from diagnostic.  Allocates 
#		relevant atomic data structure if needed.  

real procedure zcalc_dens (atl, atom, ion, assumed, diag_ratio)

#  Arguments:
pointer	atl			# IO: atomic data object list
int	atom, ion		# I: atom, ion of diagnostic
real	assumed			# I: assumed T_e
real	diag_ratio		# I: N_e sensitive diagnostic ratio
real	density			# O: computed density

#  Declarations:
pointer	at			# atomic data object
pointer	atl_get_at()		# fetch atomic data object from list
int	calc_density()		# Calculate N_e, from T_e & diagnostic
char	diag_expr[SZ_LINE]	# expression for diagnostic line ratio
int	i			# generic

errchk	at_alloc, atl_get_at, calc_density 

begin
	# Ensure that diagnostic ratio and atomic data structure are 
	# available.   
	if ( IS_INDEFR(diag_ratio) )
	    return (INDEFR)

	at = atl_get_at (atl, atom, ion)
	if (at == NULL) {
	    call at_alloc (atom, ion, at)
	    call atl_add (atl, at)
	}

	# Get expression for pre-defined diagnostic ratio.
	call set_diag_expr (at, DENSITY, diag_expr, SZ_LINE)

	i = calc_density (at, assumed, density, 1, diag_ratio, diag_expr, 
				TOLER, VERB)

	return (density)
end


#-------------------------------------------------------------------------------
#  ZCALC_TEMP -	Calculates pre-defined T_e from diagnostic.  Allocates 
#		relevant atomic data structure if needed.  

real procedure zcalc_temp (atl, atom, ion, assumed, diag_ratio)

#  Arguments:
pointer	atl			# IO: atomic data object list
int	atom, ion		# I: atom, ion of diagnostic
real	assumed			# I: assumed N_e 
real	diag_ratio		# I: T_e sensitive diagnostic ratio
real	temperature		# O: computed temperature

#  Declarations:
pointer	at			# I: atomic data object
pointer	atl_get_at()		# fetch atomic data object from list
int	calc_temperature()	# Calculate T_e, from N_e & diagnostic
char	diag_expr[SZ_LINE]	# expression for diagnostic line ratio
int	i			# generic

errchk	at_alloc, atl_get_at, calc_temperature 

begin
	# Ensure that diagnostic ratio and atomic data structure are 
	# available.   
	if ( IS_INDEFR(diag_ratio) )
	    return (INDEFR)

	at = atl_get_at (atl, atom, ion)
	if (at == NULL) {
	    call at_alloc (atom, ion, at)
	    call atl_add (atl, at)
	}

	# Get expression for pre-defined diagnostic ratio.
	call set_diag_expr (at, TEMPERATURE, diag_expr, SZ_LINE)

	i = calc_temperature (at, assumed, temperature, 1, diag_ratio, 
				diag_expr, TOLER, VERB)

	return (temperature)
end



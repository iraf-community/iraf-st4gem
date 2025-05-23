# Copyright restrictions apply - see st4gem$copyright.stsdas

include	"../atom.h"

#--------------------------------------------------------------------10 Jul 97--
.help calc_nt.x Mar97 nebular/temden
.ih
NAME
.nf
    calc_density - Calculate N_e, from T_e and diagnostic line ratio 
calc_temperature - Calculate T_e, from N_e and diagnostic line ratio 
  calc_intensity - Calculate intensity given T_e, N_e arrays & diagn expression
.fi
.ih
DESCRIPTION
Routines to calculate electron densities for an array of assumed 
temperatures, or vice versa, for a given diagnostic line ratio.  
.endhelp
#-------------------------------------------------------------------------------
#  CALC_DENSITY - 	Calculate N_e, from T_e array & diagnostic line ratio. 
#			Returns number of non-INDEF values in array.  

int procedure calc_density (at, temperature, density, n_pts, diag_ratio, 
				diag_expr, toler, verb)

# Arguments:
pointer	at			# I: atomic data structure
real	temperature[ARB]	# I: assumed electron temperature array
real	density[ARB]		# O: calculated electron density array
int	n_pts			# I: portion of density array to compute
real	diag_ratio		# I: density-sensitive diagnostic ratio
char	diag_expr[ARB]		# I: expression for density-sensitive diag ratio
real	toler			# I: fract. tolerance of temden calculation
bool	verb			# I: print diagnostic output for each iteration?
int	good_pts		# O: no. successful calculations

# Declarations:
int	i			# generic
int	status			# error code from temden routine
int	temden()		# core routine to calculate temp/density

errchk	temden

begin
	good_pts = 0
	if (IS_INDEF (diag_ratio) ) 
	    return (good_pts)

	do i = 1, n_pts {
	    status = temden (at, density[i], temperature[i], DENSITY, 
				diag_ratio, diag_expr, toler, verb)
	    if (status != OK)
		density[i] = INDEFR
	    else
		good_pts = good_pts + 1
	} 

	return (good_pts)
end


#-------------------------------------------------------------------------------
#  CALC_TEMPERATURE -	Calculate T_e, from N_e array &  diagnostic line ratio. 
#			Returns number of non-INDEF values in array.  

int procedure calc_temperature (at, density, temperature, n_pts, diag_ratio, 
				diag_expr, toler, verb)

# Arguments:
pointer	at			# I: atomic data structure
real	density[ARB]		# I: assumed electron density array
real	temperature[ARB]	# O: calculated electron temperature array
int	n_pts			# I: no. temperature values to compute
real	diag_ratio		# I: density-sensitive diagnostic ratio
char	diag_expr[ARB]		# I: expression for density-sensitive diag ratio
real	toler			# I: fract. tolerance of temden calculation
bool	verb			# I: print diagnostic output for each iteration?
int	good_pts		# O: no. successful calculations

# Declarations:
int	i			# generic
int	status			# error code from FIVEL routines
int	temden()		# core routine to calculate temp/density

errchk	temden

begin
	good_pts = 0
	if (IS_INDEF (diag_ratio) ) 
	    return (good_pts)

	do i = 1, n_pts {
	    status = temden (at, density[i], temperature[i], TEMPERATURE, 
				diag_ratio, diag_expr, toler, verb)
	    if (status != OK)
		temperature[i] = INDEFR
	    else
		good_pts = good_pts + 1
	} 

	return (good_pts)
end


#-------------------------------------------------------------------------------
#  CALC_INTENSITY -	Calculate intensity given T_e, N_e arrays. 
#			Returns number of non-INDEF values in array.  

int procedure calc_intensity (at, density, temperature, intensity, n_pts, 
				diag_expr)

# Arguments:
pointer	at			# I: atomic data structure
real	density[ARB]		# I: assumed electron density array
real	temperature[ARB]	# I: assumed electron temperature array
real	intensity[ARB]		# O: computed intensity ratio
int	n_pts			# I: no. intensity values to compute
char	diag_expr[ARB]		# I: expression for diag ratio
int	good_pts		# O: no. successful calculations

# Declarations:
int	i			# generic
real	j_ratio()		# compute specified ratio from expression

begin
	good_pts = 0
	do i = 1, n_pts {
	    iferr (call solve (at, density[i], temperature[i])) {
		intensity[i] = INDEFR
		next
	    }

	    intensity[i] = j_ratio (at, diag_expr)
	    good_pts = good_pts + 1
	} 

	return (good_pts)
end



# Copyright restrictions apply - see st4gem$copyright.stsdas

include "../at.h"
include	"../flux.h"

define	N_AVG		4	# Max abundance values to average
define	DEBUG		false

#--------------------------------------------------------------------9 Jul 97---
.help abund.x Jul97 nebular/abund
.ih
NAME
.nf
   abund - Calculates ionic abundances of ions in all I.P. zones
 ab_calc - Calculate ionic abundance from specified emission line
.fi
.endhelp
#-------------------------------------------------------------------------------
#  ABUND -	Calculates ionic abundances of all ions, relative to H+.

procedure abund2 (lfl, T_e, N_e)

#  Calling arguments
pointer	lfl			# flux list data structure
real	T_e[ARB]		# electron temperature per zone
real	N_e[ARB]		# electron density per zone

#  Declarations
real	ab_calc()		# calculate ionic abundance from em. line
pointer	at			# atomic data structure
pointer	lf			# flux data structure
int	i, j			# generic
int	n_ions			# number of ions for which to compute abundance
int	n_lines			# number of emission lines for abundance avg
int	sz_arr			# current size of val, wt arrays
pointer	wt			# relative weights
real	wt_avg()		# Wt avg of array vals, excluding INDEFs
pointer	val			# abundance from each emission line
int	zone			# ionization zone

#  Memory management
define	Val	Memr[val+$1-1]
define	Wt	Memr[wt+$1-1]

errchk	ab_calc

begin
	sz_arr = N_AVG
	call malloc (val, sz_arr, TY_REAL)
	call malloc (wt,  sz_arr, TY_REAL)

	# Operate on each ion in the flux list.
	n_ions = LFL_N(lfl)
	do i = 1, n_ions {
	    lf = LFL_A(lfl, i)
	    call at_alloc (LF_ATOM(lf), LF_ION(lf), at)
	    zone = LF_ZONE(lf)

	    n_lines = LF_N(lf)
	    if (n_lines > sz_arr) {
		call realloc (val, N_AVG, TY_REAL)
		call realloc (wt, N_AVG, TY_REAL)
		sz_arr = n_lines
	    }

	    # Calculate abundance for each line
	    do j = 1, n_lines {
		Val(j) = ab_calc (at, LF_FLUX(lf,j), LF_WAVE(lf,j), 
				N_e[zone], T_e[zone], LF_WIDTH(lf,j))
		Wt(j) = LF_WT(lf,j)
	    }

	    # Adopted abundance is weighted average of all lines.
	    LF_ABUND(lf) = wt_avg (Val(1), Wt(1), n_lines)
	}

	call at_free (at)
	call mfree (val, TY_REAL)
	call mfree (wt, TY_REAL)

	if (DEBUG)
	    call lfl_debug (lfl)
end


#-------------------------------------------------------------------------------
#  AB_CALC -	Calculate ionic abundance for specified line. 

real procedure ab_calc (at, flux, wave, dens, temp, width)

#  Arguments:
pointer	at		# I: atomic data structure
real	flux		# I: diagnostic line ratio
real	wave		# I: wavelength of emission line
real	dens		# I: electron density
real	temp		# I: electron temperature 
real	width		# I: wavelength interval
real	result		# O: result of calculation

#  Declarations:
real	emisv()		# compute emissivity from one or more lines
real	hb_emiss()	# compute H-beta emissivity

begin
	# Protect against invalid input.
	if (IS_INDEFR (flux) || IS_INDEFR (wave)) {
	    result = INDEFR
	    return (result)
	} 

	# Solve for the level populations
	iferr (call solve (at, dens, temp) ) {
	    return (INDEFR)
	}

	# Compute emissivity for all wavelengths within specified tolerance. 
	result = emisv (flux, wave, width, dens, E_TRANS(at), EMISS(at), 
			AT_NLVL(at))

	# Derive relative ionic abundance if "wave" and "flux" specified. 
	if (result > 0. && !IS_INDEF(result)) 
	    result = flux * hb_emiss(temp) / result

	else
	    result = INDEFR

	return (result)
end



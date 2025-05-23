# Copyright restrictions apply - see st4gem$copyright.stsdas

#--------------------------------------------------------------------3 Feb 97---
.help emisiv.x Feb97 nebular/fivel
.ih
NAME
.nf
    emisv - Accumulate the emissivities for specified range of lines
 hb_emiss - Compute H-beta emissivity for specified temperature.  
.fi
.endhelp
#-------------------------------------------------------------------------------
#  EMISV -	Accumulate the emissivities for specified range of lines.  

real procedure emisv (flux, wave, width, dens, energy, jl, n_lvls)

#  Arguments:
real	flux		# I: diagnostic line ratio
real	wave		# I: wavelength of emission line
real	width		# I: wavelength half-interval
real	dens		# I: input electron density
double	energy[ARB]	# I: energy above G.S. for each level
double	jl[n_lvls,ARB]	# I: volume emissivity matrix
int	n_lvls		# I: no. atomic energy levels
real	result		# O: result of calculation

#  Local variables:
real	emis		# emissivities of lines  in wavelength interval
int	i, j		# loop indexes
real	wdiff		# difference betw/requested & calcualted Waves

include	<mach.h>

define	TINY	log10 (1./MAX_REAL)
errchk	solve

begin
	if (IS_INDEFR (flux) || IS_INDEFR (wave)) {
	    result = INDEFR
	    return (result)
	}

	# Derive wavelengths & emissivities from all nearby lines. 
	emis = 0.
	do i = 2, n_lvls {
	    do j = 1, i-1 {
	    	wdiff = abs (wave - 1.D+0 / (energy[i] - energy[j]) )
	    	if (wdiff <= width && !IS_INDEF(jl[i,j]) ) {
		    emis = emis + jl[i,j]
		}
	    }
	}

	# Derive relative ionic abundance if "wave" and "flux" specified.  
	# Protect against underflow. 
	if (emis > 0.) {
	    if (log10(emis) - log10(dens) > TINY)
	    	result = emis / dens
	    else
	    	result = INDEFR
	} else 
	    result = 0.

	return (result)
end


#-------------------------------------------------------------------------------
#  HB_EMISS -	Compute H-beta emissivity for specified temperature.  

real procedure hb_emiss (temp)

#  Arguments:
real	temp		# I: electron temperature 
real	j_hb		# O: result of calculation

#  Declarations:
real	t4		# temperature / 10,000

begin
	if (temp < 0. || temp > 1.e+6) 
	    return (INDEFR)

	# Use Brocklehurst's interpolation formula for H-beta emissivity. 
	# Accurate within 4% for density <= 1.e+6/cm^3. 
	t4   = temp * 1.e-4
	j_hb = 1.387 / (t4 ** 0.983) / 10. ** (0.0424 / t4) * 1.e-25

	return (j_hb)
end



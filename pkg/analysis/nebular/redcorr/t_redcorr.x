# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<error.h>
include	<mach.h>
include	"../fivel.h"
include	"../neberr.h"

#--------------------------------------------------------------------19 Feb 97--
.help t_redcorr.x Mar95 nebular
.ih
NAME
.nf
t_redcorr - Pperform interstellar redding correction for a specified line
    dered - Correct input fluxes for interstellar extinction. 
     isef - Compute the interstellar extinction function 
.fi
.ih
DESCRIPTION
.endhelp
#-------------------------------------------------------------------------------
#  T_REDCORR -	Perform interstellar redding correction for a specified line.  

procedure t_redcorr ()

#  Local variables
char	keywd[SZ_FNAME]		# matched keyword
real	c_ext			# log extinction @ H-beta
real	flux			# line flux relative to H-beta
int	rfunc			# redding function code from "fivel.h"
real	val			# input electron temperature or density
real	wave			# wavelength of line (Angstroms)

#  Functions used:
real	clgetr()		# fetch CL parameter TY_REAL
int	clgwrd()		# fetch a keyword from an enumerated string

errchk	dered

begin

	# Get the task parameters
	wave  = clgetr ("wave")
	flux  = clgetr ("flux")
	c_ext = clgetr ("c_ext")

	rfunc = clgwrd ("red_func", keywd, SZ_FNAME, EXTN_MODEL)
	if (rfunc <= 0) 
	    call error (BAD_EXTN_LAW, "Invalid reddening function")
	call clpstr ("red_func", keywd)

	# Assure valid wavelength 
	if (IS_INDEFR(wave)) 
	    call error (BAD_EXTN_WAVE, "Wavelength must exceed 1000 Ang")

	# Perform calculation & write out result
	call dered (flux, wave, val, 1, rfunc, c_ext)

	call clputr ("result", val)

	# Write results to STDOUT
	call strupr (keywd) 
	call printf (
	"# Reddening correction using %s function:\n  Flux: %7g\n")
	    call pargstr (keywd)
	    call pargr (val)

end


#-------------------------------------------------------------------------------
#  DERED - 	Correct an array of input fluxes for interstellar extinction. 
#		Checks for physical flux values & sets to INDEF if flux is 
#		INDEF. 

procedure dered (flux, wave, val, npts, red_law, c)

#  Arguments:
real	flux[ARB]	# I: array of fluxes to be corrected
real	wave[ARB]	# I: wavelength array, Angstroms
real	val[ARB]	# O: corrected fluxes
int	npts		# I: size of arrays
int	red_law		# I: choice of redding law from "fivel.h"
real	c		# I: logarithmic extinction at H-beta (4861 Ang)

#  Local variables:
int	i		# generic
real	r_v		# ratio of total to selective extinction (@V)

errchk	isef

begin

	# Find extinction correction, assuming R_v = 3.1
	r_v = 3.1

	# Assure valid fluxes.
	do i = 1, npts {
	    if ( flux[i] < 0. || IS_INDEFR(flux[i]) )
	    	call error (1, "Invalid flux for extinction correction")

	}

	call isef (wave, val, npts, red_law, c, r_v)
	call amulr (val, flux, val, npts)

end


#-------------------------------------------------------------------------------
#  ISEF -	Compute the interstellar extinction function and renormalize 
#		to the relative extinction at H-beta (4861 Ang).

procedure isef (wave, val, npts, red_law, c, r_v)

#  Calling arguments:
real	wave[ARB]	# I: wavelength array, Angstroms
real	val[ARB]	# O: computed correction array
int	npts		# I: size of arrays
int	red_law		# I: code for I.S. redding law 
real	c		# I: logarithmic extinction at H-beta (4861 Ang)
real	r_v		# I: ratio of total to selective extinction (@V)

#  Local variables:
int	i		# generic

begin
	if ( IS_INDEFR (c) ) 
	    call error (4, "Invalid extinction constant") 

	switch (red_law) {
	case GAL:
	    call gal_redlaw (wave, val, npts)
	    call adivkr (val, 1.169930, val, npts)
	    call asubkr (val, 1.0,      val, npts)

	case CCM:
	    if ( IS_INDEFR (r_v) ) 
	    	call error (4, "Bad inputs for extinction correction") 
	    call ccm_redlaw (wave, val, npts, r_v)
	    call adivkr (val, 1.164172, val, npts)
	    call asubkr (val, 1.0,      val, npts)

	case LMC:
	    call lmc_redlaw (wave, val, npts)
	    call adivkr (val, 1.157702, val, npts)
	    call asubkr (val, 1.0,      val, npts)

	case SMC:
	    call smc_redlaw (wave, val, npts)
	    call adivkr (val, 1.144282, val, npts)
	    call asubkr (val, 1.0,      val, npts)

	case JBK:
	    call jbk_redlaw (wave, val, npts)

	default:
	    call error (4, "Invalid extinction law") 

	}

	do i = 1, npts 
	    val[i] = 10. ** (c * val[i])
end



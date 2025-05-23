# Copyright restrictions apply - see st4gem$copyright.stsdas

include	"../at.h"
include	"../atom.h"
include	"../neberr.h"

define	DEBUG	false

#---------------------------------------------------------------------6 Feb 97--
.help t_ionic.x May96 nebular/fivel
.ih
NAME
      t_ionic - Task to compute line emissivities & ionic abundances. 
.endhelp
#-------------------------------------------------------------------------------
#  T_IONIC -	Task to performs ionic abundance calculations, given nebular 
#		density, temperature, and observed line wavelength & flux.  
#		Note: flux must be in units of F(H-beta) = 100. 

procedure t_ionic ()

#  Local variables
pointer	at			# atomic data structure
int	atno			# atomic number
real	dens			# input electron density
char	errmsg[SZ_LINE]		# text of error message
real	flux			# input line ratio
int	i			# generic
char	keywd[SZ_FNAME]		# matched keyword
real	result			# result of calculation
int	stage			# ionization stage--i.e., spectrum number
real	temp			# input electron temperature 
bool	verb			# print level populations & critical densities?
real	wave			# wavelength of emission line
real	width			# wavelength half-interval for emission lines

#  Functions used:
bool	clgetb()		# fetch CL parameter TY_BOOL
int	clgeti()		# fetch CL parameter TY_INT
real	clgetr()		# fetch CL parameter TY_REAL
int	clgwrd()		# fetch a keyword from an enumerated string
int	get_atom_typ()		# get atomic number from input string

errchk	abund_out, pari, parii, solve

begin

	# Get the task parameters; update the "temden", "atom", and 
	# "ion_stage" task parameters
	temp  = clgetr ("temperature")
	dens  = clgetr ("density")
	wave  = clgetr ("wave")
	width = clgetr ("wv_toler")
	flux  = clgetr ("flxratio") 
	if ( !IS_INDEFR (flux) )
	    flux = flux / 100.
	verb = clgetb ("verbose")

	i = clgwrd ("atom", keywd, SZ_FNAME, ATOM)
	call clpstr ("atom", keywd)
	atno = get_atom_typ (keywd)

	stage = clgeti ("spectrum")
	call clputi ("spectrum", stage)

	# Assure valid input temp/dens
	if (IS_INDEFR(temp)) 
	    call error (BAD_T_E, "Must provide electron temperature!")
	else
	    call clputr ("temperature", temp)

	if (IS_INDEFR(dens)) 
	    call error (BAD_N_E, "Must provide electron density!")
	else
	    call clputr ("density", dens)

	# Allocate atomic data structure.  
	call at_alloc (atno, stage-1, at)

	# Check for in-bounds T_e. 
	if ( temp < AT_TE_MIN(at) || temp > AT_TE_MAX(at) ) {
	    call clputr ("result", INDEFR)
	    call sprintf (errmsg, SZ_LINE, 
		"T_e out of bounds for this ion. Use %8.1f <= T_e <= %8.1f")
	    call pargr (AT_TE_MIN(at))
	    call pargr (AT_TE_MAX(at))
	    call at_free (at)
	    call error (ABUND_FAILED, errmsg)
	} 

	iferr (call solve (at, dens, temp)) {
	    call clputr ("result", INDEFR)
	    call at_free (at)
	    call sprintf (errmsg, SZ_LINE, "Emissivity calculation failed")
	    call error (ABUND_FAILED, errmsg)
	} 

	if (DEBUG)
	    call at_debug (at)

	# Write out results.
	call ion_out (at, flux, wave, dens, temp, result, width, verb)

	call at_free (at)
end



# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<evexpr.h>
include	<error.h>
include	"../atom.h"
include	"../at.h"
include	"../neberr.h"

define	DEBUG	false

#--------------------------------------------------------------------10 Jul 97--
.help t_temden.x May96 nebular/temden
.ih
NAME
.nf
  t_temden -- SPP driver for temperature/density calculation 
    out_td -- Write TEMDEN output to STDOUT 
.fi
.endhelp
#-------------------------------------------------------------------------------
#  T_TEMDEN -	SPP driver for five-level atom program.  Performs nebular 
#		density and temperature calculations.  

procedure t_temden ()

#  Local variables
pointer	at			# atomic data structure
int	atom			# atomic number
char	diag_expr[SZ_LINE]	# expression for density-sensitive diag ratio
int	diag_type		# diagnostic type: "density|temperature"
pointer	ep			# pointer to expression structure
char	expstr[SZ_FNAME]	# string with expression for line ratio 
int	i			# generic
int	ion			# ion--i.e., spectrum number
char	keywd[SZ_FNAME]		# matched keyword
real	ratio			# diagnostic line ratio
real	result			# result of calculation
real	toler			# fract. tolerance on result
char	transition[SZ_FNAME]	# expression for transition of interest
real	assume			# input electron temperature or density
bool	verb			# print diagnostic output for each iteration?

#  Functions used:
int	calc_density()		# calculate N_e
int	calc_temperature()	# calculate T_e
bool	clgetb()		# fetch CL parameter TY_BOOL
int	clgeti()		# fetch CL parameter TY_INT
real	clgetr()		# fetch CL parameter TY_REAL
int	clgwrd()		# fetch a keyword from an enumerated string
pointer	evexpr()		# evaluate an algebraic expression
int	get_atom_typ()		# get atomic number from string variable
bool	streq()			# are two strings equal?

errchk	calc_density, calc_temperature, evexpr, out_td

begin
	# Get the task parameters; update the "option", "atom", and "ion_stage" 
	# task parameters
	diag_type = clgwrd ("option", keywd, SZ_FNAME, TEMDEN)
	call clpstr ("option", keywd)

	# Evaluate expression for line flux ratio
	call clgstr ("flxratio", expstr, SZ_FNAME)
	ep = evexpr (expstr, NULL, NULL)
	switch (O_TYPE(ep)) {
	case TY_REAL:
	    ratio = O_VALR(ep)
	default:
	    call error (BAD_RATIO, 
		"Expression must involve floating-point numbers")
	}
	call mfree (ep, TY_STRUCT)

	i = clgwrd ("atom", keywd, SZ_FNAME, ATOM)
	call clpstr ("atom", keywd)
	atom = get_atom_typ (keywd)

	# Note: ion number is zero-indexed, but task parameter is one-indexed.
	ion = clgeti ("spectrum")
	call clputi ("spectrum", ion)
	ion = ion - 1

	assume = clgetr ("assume")
	call clputr ("assume", assume)

	verb = clgetb ("verbose")

	# Assure valid line ratio and input temperature/density.
	if (IS_INDEFR(ratio)) 
	    call error (BAD_RATIO, "Must provide diagnostic line ratio")

	if (IS_INDEFR(assume)) 
	    call error (BAD_NT, "Must provide electron temperature or density")

	# Allocate/initialize the atomic data structure. 
	call at_alloc (atom, ion, at)
	if (DEBUG)
	    call at_debug (at)

	# Fetch transition description.
	call clgstr ("transition", transition, SZ_FNAME)
	call strupr (transition)
	if ( streq (transition, "DEFAULT") || streq (transition, EOS) ) 
	    call set_diag_expr (at, diag_type, diag_expr, SZ_FNAME)

	else
	    call strcpy (transition, diag_expr, SZ_FNAME)

	# Compute temperature or density.
	toler = 1.e-5
	if (diag_type == DENSITY) {
	    i = calc_density (at, assume, result, 1, ratio, diag_expr, 
				toler, verb)
	    if (IS_INDEFR (result)) 
	    	call printf ("Density calculation failed\n")

	} else if (diag_type == TEMPERATURE) {
	    i = calc_temperature (at, assume, result, 1, ratio, diag_expr, 
				toler, verb)
	    if (IS_INDEFR (result)) 
	    	call printf ("Temperature calculation failed\n")
	}
	if (DEBUG)
	    call at_debug (at)

	# Write results to pset and STDOUT. 
	call clputr ("result", result)
	if ( !IS_INDEFR (result) )
	    call out_td (at, ratio, result, diag_type, transition)

	call at_free (at)
end


#-------------------------------------------------------------------------------
#  OUT_TD -	Write TEMDEN output to STDOUT. 

procedure out_td (at, ratio, result, diag_type, transition)

#  Calling arguments:
pointer	at			# atomic data structure
real	ratio			# diagnostic line ratio
real	result			# result of calculation
int	diag_type		# diagnostic type: "density|temperature"
char	transition[ARB]		# transition description

#  Declarations:
int	get_atom_str()		# Get the atom/ion string
int	get_ion_str()		# Get the atomic transition description
char	ion_string[SZ_FNAME]	# description of atomic transition
int	nchars			# length of ion_string
bool	streq()			# are two strings equal?
char	sx[SZ_LINE]		# generic

begin
	# Fetch the appropriate ion string. 
	nchars = get_ion_str (at, diag_type, ion_string, SZ_FNAME)
	if ( !(streq (transition, "DEFAULT") || streq (transition, EOS)) ) {
	    nchars = get_atom_str (AT_ATOM(at), AT_ION(at), sx, SZ_LINE)
	    call sprintf (ion_string, SZ_LINE, "%s: %s")
		call pargstr (sx)
		call pargstr (transition)
	}

	# Format and print the transition description...
	if (nchars > 0) {
	    call printf ("%s %s = %6g\n")
	    if (diag_type == DENSITY)
	    	call pargstr ("Density ratio ")

	    else if (diag_type == TEMPERATURE)
	    	call pargstr ("Temperature ratio ")

	    call pargstr (ion_string)
	    call pargr (ratio)

	# ...unless the specified transition is invalid.
	} else {
	    call printf ("No nebular diagnostics available for that ion\n")
		return
	}
	call flush (STDOUT)

	# Format & print the calculated ratio.
	if (diag_type == DENSITY) 
	    call printf ("Density = %7g /cm^3\n")

	else if (diag_type == TEMPERATURE) 
	    call printf ("Temperature = %6.1f K\n")
	
	call pargr (result)
end

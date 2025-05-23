# Copyright restrictions apply - see st4gem$copyright.stsdas

include	"../atom.h"

#---------------------------------------------------------------------6 Feb 98--
.help valid.x Aug96 nebular
.ih
NAME
   diag_valid - Ensure a supported diagnostic
.endhelp
#-------------------------------------------------------------------------------
#  DIAG_VALID - Is the requested diagnostic supported for the atom/ion?  
#		Diagnostic types are: emissivity (0), density (1), and 
#		temperature (2). 

bool procedure diag_valid (at, diag_type)

#  Arguments:
pointer	at		# I: atomic data structure
bool	valid		# O: is ionum valid?
int	diag_type	# I: diagnostic type: emissivity|density|temperature 

#  Local variables:
int	atom		# atomic number
int	spect		# spectrum number

define	EMISSIVITY	0

include	"../at.h"

begin
	atom  = AT_ATOM(at)
	spect = AT_ION(at)

	valid = false

	switch (atom) {
	case CARBON:
	    if ( spect == 0 && (diag_type != DENSITY ) )
	    	valid = true
	    else if ( spect == 1 && (diag_type != TEMPERATURE ) )
	    	valid = true
	    else if ( spect == 2 && (diag_type != TEMPERATURE) )
	    	valid = true

	case NITROGEN:
	    if (spect == 0)
		valid = true
	    else if ( spect == 1 && (diag_type != DENSITY ) )
		valid = true
	    else if ( spect == 2 && (diag_type != TEMPERATURE) )
		valid = true
	    else if ( spect == 3 && (diag_type != TEMPERATURE) )
		valid = true

	case OXYGEN:
	    if (spect == 0 && (diag_type != DENSITY) ) 
	    	valid = true
	    else if (spect == 1) 
	    	valid = true
	    else if (spect == 2 && (diag_type != DENSITY) ) 
	    	valid = true
	    else if (spect == 3 && (diag_type != TEMPERATURE) ) 
	    	valid = true
	    else if (spect == 4 && (diag_type != TEMPERATURE) ) 
	    	valid = true

	case NEON:
	    if (spect == 2 && (diag_type != DENSITY) ) 
	    	valid = true
	    else if (spect == 3) 
	    	valid = true
	    else if (spect == 4 && (diag_type != DENSITY) ) 
	    	valid = true
	    else if (spect == 5 && (diag_type != TEMPERATURE) ) 
	    	valid = true

	case SODIUM:
	    if (spect == 3 && (diag_type != DENSITY) ) 
	    	valid = true
	    else if (spect == 5 && (diag_type != DENSITY) ) 
	    	valid = true

	case MAGNESIUM:
	    if (spect == 4 && (diag_type != DENSITY) ) 
	    	valid = true

	    if (spect == 6 && (diag_type != DENSITY) ) 
	    	valid = true

	case ALUMINUM:
	    if (spect == 1)
	    	valid = true

	case SILICON:
	    if (spect == 1 && (diag_type != TEMPERATURE) ) 
	    	valid = true
	    else if (spect == 2) 
	    	valid = true

	case SULFUR:
	    if (spect == 1) 
	    	valid = true
	    else if (spect == 2 && (diag_type != DENSITY) ) 
	    	valid = true
	    else if (spect == 3 && (diag_type != TEMPERATURE) ) 
	    	valid = true

	case CHLORINE:
#	    if (spect == 1 && (diag_type != DENSITY) )
#	    	valid = true
	    if (spect == 2) 
	    	valid = true
	    else if (spect == 3 && (diag_type != DENSITY) ) 
	    	valid = true

	case ARGON:
	    if (spect == 2 && (diag_type != DENSITY) ) 
	    	valid = true
	    else if (spect == 3) 
	    	valid = true
	    else if (spect == 4 && (diag_type != DENSITY) ) 
	    	valid = true

	case POTASSIUM:
	    if (spect == 3 && (diag_type != DENSITY) )
	    	valid = true
	    else if (spect == 4) 
	    	valid = true

	case CALCIUM:
	    ;
#	    if (spect == 4 && (diag_type != DENSITY) ) 
	}

	return (valid)
end



# Copyright restrictions apply - see st4gem$copyright.stsdas

include	"../at.h"
include	"../atom.h"

#---------------------------------------------------------------------6 Feb 97--
.help diagnostic.x Jan96 nebular/fivel
.ih
NAME
      j_ratio - Calculate the specified ratio from the line emissivity matrix
set_diag_expr - Set the diagnostic expression for pre-defined ratios
.endhelp
#-------------------------------------------------------------------------------
#  J_RATIO --	Calculate the specified ratio from the line emissivity matrix.
#		For now, support diagnostics for a single atom at a time.  

real procedure j_ratio (at, diag_expr)

#  Arguments:
pointer	at			# I: atomic data structure
char	diag_expr[ARB]		# I: expression for diagnostic ratio
real	rcal			# O: line ratio

#  Functions used:
real	at_evdiag()

begin
	rcal = at_evdiag (at, diag_expr)
	return (rcal)
end


#-------------------------------------------------------------------------------
#  SET_DIAG_EXPR - Set the diagnostic expression for pre-defined ratios.

procedure set_diag_expr (at, diag, expression, max_char)

#  Arguments:
pointer	at			# I: atomic data structure
int	diag			# I: diagnostic type (density|temperature)
char	expression[ARB]		# O: expression for diagnostic ratio
int	max_char		# I: max size of expression string

#  Declarations:
int	ion			# ion

begin
	ion = AT_ION(at)

	# Most diagnostic ratios are set by the ground-state electron 
	# configurations. 
	switch (AT_GSCONFIG(at)) {
	case GS_S2:
	    if (diag == DENSITY)
	    	call strcpy ("J(4,1) / J(3,1)", expression, max_char)

	    else if (diag == TEMPERATURE) 
	    	call strcpy ("(J(4,1) + J(3,1)) / J(5,1)", expression, max_char)

	case GS_P1:
	    if (diag == DENSITY)
	    	call strcpy ("J(5,2) / J(4,2)", expression, max_char)

	case GS_P2:
	    if (diag == TEMPERATURE)
	    	call strcpy ("(J(4,2) + J(4,3)) / J(5,4)", expression, max_char)

	case GS_P3:
	    if (diag == DENSITY)
	    	call strcpy ("J(3,1) / J(2,1)", expression, max_char)

	    else if (diag == TEMPERATURE) 
	    	    call strcpy ("(J(2,1) + J(3,1)) / (J(4,1) + J(5,1))", 
				expression, max_char)

	case GS_P4:
	    if (diag == TEMPERATURE)
	    	call strcpy ("(J(4,1) + J(4,2)) / J(5,4)", expression, max_char)

	default:
	    call strcpy (" ", expression, max_char)
	}

	# Now handle the exceptions for p^3 ions.
	switch (AT_ATOM(at)) {
	case NITROGEN:
	    if (ion == 0) {
		if (diag == TEMPERATURE) 
	    	call strcpy ("(J(2,1)+J(3,1)) / (J(5,3)+J(5,2)+J(4,3)+J(4,2))", 
				expression, max_char)
	    }

	case OXYGEN:
	    if (ion == 1) {
		if (diag == TEMPERATURE) 
	    	call strcpy ("(J(2,1)+J(3,1)) / (J(5,3)+J(5,2)+J(4,3)+J(4,2))", 
				expression, max_char)
	    }
	}
end



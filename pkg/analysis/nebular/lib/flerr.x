# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<error.h>
include	"../neberr.h"

#--------------------------------------------------------------------6 Feb 97---
.help flerr.x Mar96 nebular/fivel
.ih
NAME
  fl_err - Write FIVEL error status info to STDOUT. 
.endhelp
#-------------------------------------------------------------------------------
#  FL_ERR -	Write FIVEL error status info to STDOUT. 

procedure fl_err (status)

#  Calling argument:
int	status			# return status of FIVEL routine

include	"../flerr.com"

begin
	switch (status) {
	case NOT_CONVERGENT:
	    call printf ("WARNING: Did not converge after 100 iterations.\n")

	case ARITHMETIC_ERROR:
	    call printf ("Divide by zero in routine SOLVE.\n")

	case LR_OUT_OF_BOUNDS:
	    call printf ("EARLY EXIT -- Line ratio out of bounds.\n")
	    call printf (
		"The entered RATIO: %10.5e predicts unreasonable conditions.\n")
		    call pargr (val1)
	    call printf (
		"The allowable RATIO range for atom %2d is: %9.5e -- %9.5e\n")
		    call pargi (int(val2))
		    call pargr (val3)
		    call pargr (val4)

	case LR_NEAR_LIMIT:
	    call printf ("WARNING: Line ratio near useful limit.\n")

	case ION_INVALID:
	    call printf ("Specified atom/spectrum diagnostic unavailable.\n")
	}

end



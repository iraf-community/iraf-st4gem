include "nlfit.h"

# NL_MINIMIZE -- Call the appropriate minimization algorithm.

int procedure nl_minimize (nl)

pointer nl		# i: Curve descriptor.

#--
int	nl_stati(), nl_amoeba(), nl_marq()

errchk	nl_amoeba, nl_marq

begin
	switch (nl_stati (nl, "algorithm")) {

	case AMOEBA:
	    return (nl_amoeba (nl))

	case MARQUARDT:
	    return (nl_marq (nl))

	default:
	    call error (0, "Non-supported minimization method.\n")
	}
end

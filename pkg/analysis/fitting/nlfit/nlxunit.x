include "nlfit.h"

define	LIMMETER     5.e-6	# X lower than this implies wavlength in meter.
define	LIMCM        1.e-3	# Limit for cm.
define	LIMANG1	      300.	# Limits for Angstrom physical units.
define	LIMANG2	    50000.	# If x greater than this, then Hz.

# NL_XUNIT -- Chooses a physical unit for Planck function evaluation.

int procedure nl_xunit (x)

real	x

#--
begin
	if (x < LIMMETER)
	    return (METER)
	else if ((x >= LIMMETER) && (x < LIMCM))
	    return (CM)
	else if ((x >= LIMANG1) && (x < LIMANG2))
	    return (ANGSTROM)
	else if (x >= LIMANG2)
	    return (HZ)
	else {
	    call eprintf ("\nWarning: no discernible physical units. ")
	    return (ERR)
	}
end
                          
             

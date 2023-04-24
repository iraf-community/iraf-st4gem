include "addmasks.h"

# MAPTRUE -- Compute the value of true (all flags set)
#
# B.Simon	16-Jul-93	Original
# B.Simon	19-Oct-93	Modified to work around compiler bug

int procedure maptrue (nflag)

int	nflag		# i: Number of flag values
#--
int	iflag, tval

begin
	if (nflag == 0) {
	    iflag = MAX_FLAGS
	} else {
	    iflag = nflag
	}

	tval = 2 ** iflag - 1
	return (tval)
end

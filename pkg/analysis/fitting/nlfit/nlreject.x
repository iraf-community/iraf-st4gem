include "nlfit.h"

# NL_REJECT -- Procedure to eliminate a single datapoint from the data set.

procedure nl_reject (nl, i)

pointer	nl		# i: curve descriptor
int	i		# i: index of rejected point

#--
begin
	Memb[NL_REJFLAG(nl)+i-1] = true
	NL_NREJECT(nl) = NL_NREJECT(nl) + 1
end


# NL_CLEAR -- Initialize rejection, clearing all data point flags.

procedure nl_clear (nl)

pointer	nl		# i: curve descriptor

#--
int	i

begin
	NL_NREJECT(nl) = 0
	do i = 0, NL_NPTS(nl) - 1
	    Memb[NL_REJFLAG(nl)+i] = false
end

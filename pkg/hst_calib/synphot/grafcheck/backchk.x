include "../grflist.h"

# BACKCHK -- Check the instrument graph for rows with null values
#
# B.Simon	15-Aug-88	First Code

procedure backchk (gp)

pointer	gp		# i: Pointer to graph structure
#--
bool	first
int	list

string	back_intro  "The following rows have backward pointing nodes:\n\n"
string	back_fmt    "\"%s\"%25t\"%s\"%40t%10d%10d\n"

begin
	first = true

	do list = 1, GRF_SIZE(gp) {

	    if (GRF_INNODE(gp,list) >= GRF_OUTNODE(gp,list)) {

		if (first) {
		    first = false
	     	    call putline (STDOUT, back_intro)
		}

		call printf (back_fmt)
		    call pargstr (GRF_COMPNAME(gp,list))
		    call pargstr (GRF_KEYNAME(gp,list))
		    call pargi (GRF_INNODE(gp,list))
		    call pargi (GRF_OUTNODE(gp,list))
	    }
	}

	if (! first)
	    call putline (STDOUT, "\n")

end

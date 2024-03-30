include "../grflist.h"

# NULLCHK -- Check the instrument graph for rows with null values
#
# B.Simon	15-Aug-88	First Code

procedure nullchk (gp)

pointer	gp		# i: Pointer to graph structure
#--
bool	first
int	list

string	null_intro  "The following rows have undefined values:\n\n"
string	null_fmt    "\"%s\"%25t\"%s\"%40t%10d%10d\n"

begin
	first = true

	do list = 1, GRF_SIZE(gp) {

	    if (IS_INDEFI (GRF_INNODE(gp,list)) ||
		IS_INDEFI (GRF_OUTNODE(gp,list)) ||
		GRF_KEYNAME(gp,list) == EOS ||
		GRF_COMPNAME(gp,list) == EOS 	 ) {

		if (first) {
		    first = false
	     	    call putline (STDOUT, null_intro)
		}

		call printf (null_fmt)
		    call pargstr (GRF_COMPNAME(gp,list))
		    call pargstr (GRF_KEYNAME(gp,list))
		    call pargi (GRF_INNODE(gp,list))
		    call pargi (GRF_OUTNODE(gp,list))
	    }
	}

	if (! first)
	    call putline (STDOUT, "\n")

end

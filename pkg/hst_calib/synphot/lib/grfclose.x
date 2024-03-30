include	"../grflist.h"

# GRF_CLOSE -- Deallocate memory used by graph list
#
# B.Simon	12-Aug-88	First Code

procedure grf_close (gp)

pointer	gp		# io: Pointer to graph list data structure
#--

errchk	mfree

begin
	if (gp == NULL)
	    return

	call mfree (GRF_KEYARRAY(gp), TY_CHAR)
	call mfree (GRF_CMPARRAY(gp), TY_CHAR)
	call mfree (GRF_OUTARRAY(gp), TY_INT)
	call mfree (GRF_INARRAY(gp), TY_INT)

	call mfree (gp, TY_INT)
	gp = NULL

end

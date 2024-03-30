include	"../adjlist.h"

# ADJ_CLOSE -- Deallocate memory used by adjacency list
#
# B.Simon	21-Jul-88	First Code

procedure adj_close (gp)

pointer	gp		# io: Pointer to adjacency list
#--
int	list, maxlist

errchk	mfree

begin
	if (gp == NULL)
	    return

	maxlist = ADJ_SIZE(gp)

	if (ADJ_NODPTR(gp) != NULL) {
	    do list = 1, maxlist
		call mfree (ADJ_NODARY(gp,list), TY_INT)
	    call mfree (ADJ_NODPTR(gp), TY_INT)
	}

	if (ADJ_NAMPTR(gp) != NULL) {
	    do list = 1, maxlist
		call mfree (ADJ_NAMARY(gp,list), TY_CHAR)
	    call mfree (ADJ_NAMPTR(gp), TY_INT)
	}

	if (ADJ_NUMPTR(gp) != NULL)
	    call mfree (ADJ_NUMPTR(gp), TY_INT)

	call mfree (gp, TY_INT)
	gp = NULL

end

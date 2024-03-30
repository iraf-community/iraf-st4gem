include "../grflist.h"

# AMBIGCHK -- Check the instrument graph for rows with ambiguous keywords
#
# B.Simon	15-Aug-88	First Code

procedure ambigchk (gp)

pointer	gp		# i: Pointer to graph structure
#--
bool	first
int	lastnode, ilist, jlist, klist, mlist
pointer	link

string	ambig_intro  "The following rows have ambiguous keywords:\n\n"
string	ambig_fmt    "\"%s\"%25t\"%s\"%40t%10d%10d\n"

bool	streq()

begin
	# Allocate an array to hold linked lists of matching nodes

	call calloc (link, GRF_SIZE(gp), TY_INT)

	lastnode = 0
	do ilist = 1, GRF_SIZE(gp) {

	    if (GRF_INNODE(gp,ilist) != lastnode) {
		lastnode = GRF_INNODE(gp,ilist)
		klist = ilist
	    }

	    # If two rows have the same innode and the same keyword, the
	    # graph is ambiguous. Save matching rows as a linked list to
	    # be printed later.

	    do jlist = klist, ilist-1 {

		if (streq (GRF_KEYNAME(gp,jlist), GRF_KEYNAME(gp,ilist))) {

		    mlist = jlist
		    while (Memi[link+mlist-1] != 0)
			mlist = Memi[link+mlist-1]
		    Memi[link+mlist-1] = ilist

		    break
		}

	    }

	}

	first = true
	do ilist = 1, GRF_SIZE(gp) {

	    # Print ambiguous rows. Zero out the linked list so that
	    # each row is only printed once.

	    if (Memi[link+ilist-1] != 0) {

		if (first) {
		    first = false
	     	    call putline (STDOUT, ambig_intro)
		}

		call printf (ambig_fmt)
		    call pargstr (GRF_COMPNAME(gp,ilist))
		    call pargstr (GRF_KEYNAME(gp,ilist))
		    call pargi (GRF_INNODE(gp,ilist))
		    call pargi (GRF_OUTNODE(gp,ilist))

		klist = ilist
		jlist = Memi[link+ilist-1]
		Memi[link+klist-1] = 0

		while (jlist != 0) {

		    call printf (ambig_fmt)
			call pargstr (GRF_COMPNAME(gp,jlist))
			call pargstr (GRF_KEYNAME(gp,jlist))
			call pargi (GRF_INNODE(gp,jlist))
			call pargi (GRF_OUTNODE(gp,jlist))

		    klist = jlist
		    jlist = Memi[link+jlist-1]
		    Memi[link+klist-1] = 0
		}

	    	call putline (STDOUT, "\n")
	    }
	}

	call mfree (link, TY_INT)
end

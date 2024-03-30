include	"../grflist.h"

# CNECTCHK -- Check to see if graph is entirely connected
#
# B.Simon	15-Aug-88	First Code

procedure cnectchk (gp)

pointer	gp		# i: Pointer to graph structure
#--
int	list, trip, maxtrip
pointer	visit

string	cnect_intro  "The following rows are not connected to the graph:\n\n"
string	cnect_fmt    "\"%s\"%25t\"%s\"%40t%10d%10d\n"

begin
	# Perform a depth first search on the graph
	# Mark each set of connected nodes with the trip number

	call calloc (visit, GRF_SIZE(gp), TY_INT)

	trip = 1
	do list = 1, GRF_SIZE(gp) {
	    if (Memi[visit+list-1] == 0) {
		if (GRF_INNODE(gp,list) != GRF_INNODE(gp,1))
		    trip = trip + 1

		call cnect_dfs (gp, list, trip, Memi[visit])
	    }
	}
	maxtrip = trip

	# Print all nodes which were not visited in the first trip

	if (trip > 1)
	    call putline (STDOUT, cnect_intro)

	do trip = 2, maxtrip {
	    do list = 1, GRF_SIZE(gp) {

		if (Memi[visit+list-1] == trip) {
		    call printf (cnect_fmt)
			call pargstr (GRF_COMPNAME(gp,list))
			call pargstr (GRF_KEYNAME(gp,list))
			call pargi (GRF_INNODE(gp,list))
			call pargi (GRF_OUTNODE(gp,list))
		}
	    }
	    call putline (STDOUT, "\n")
	}

	call mfree (visit, TY_INT)
end

# CNECT_DFS -- Visit all connected nodes using a depth first search (dfs)

procedure cnect_dfs (gp, root, trip, visit)

pointer	gp		#  i: Pointer to graph structure
int	root		#  i: Root node in set of connected nodes
int	trip		#  i: Trip number
int	visit[ARB]	# io: Array containing trip number of node
#--
bool	found
int	top, ilist, jlist, inode, jnode
pointer	stack

bool	strne()

begin
	call malloc (stack, GRF_SIZE(gp), TY_INT)

	# Put the root row on top of the stack

	top = 0
	Memi[stack] = root
	visit[root] = trip

	while (top >= 0) {

	    ilist = Memi[stack+top]
	    inode = GRF_OUTNODE(gp,ilist)

	    # Look for a child of the row on top of the stack

	    found = false
	    if (!IS_INDEFI (inode)) {

		do jlist = ilist+1, GRF_SIZE(gp) {

		    jnode = GRF_INNODE(gp,jlist)
		    if (jnode > inode)
			break

		    if (jnode == inode && visit[jlist] == 0) {
			found = true
			top = top + 1
			Memi[stack+top] = jlist
			visit[jlist] = trip
			break
		    }

		}
	    }

	    # Find another row containing the same component as
	    # the row on top of the stack

	    if (! found) {
		top = top - 1
		do jlist = ilist+1, GRF_SIZE(gp) {

		    if (GRF_INNODE(gp,ilist) != GRF_INNODE(gp,jlist) ||
			strne (GRF_COMPNAME(gp,ilist),
			       GRF_COMPNAME(gp,jlist))			)
			break

		    visit[jlist] = trip
		    if (GRF_OUTNODE(gp,ilist) != GRF_OUTNODE(gp,jlist)) {
			top = top + 1
			Memi[stack+top] = jlist
			break
		    }
		}
	    }

	}

	call mfree (stack, TY_INT)
end

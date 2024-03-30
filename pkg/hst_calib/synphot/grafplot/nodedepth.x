include	"../adjlist.h"

# NODEDEPTH -- Calculate the depth of nodes in a graph
#
# The depth of a given node in a graph is the number of nodes that must be
# passed through from the starting node to reach the given node, inclusive
# of the given node.
#
# B.Simon	26-Jul-88	First Code

procedure nodedepth (gp, start, depth)

pointer	gp		# i: Pointer to adjacency list structure
int	start		# i: Starting node for depth calculation
real	depth[ARB]	# o: Array containing node depths
#--
int	root
pointer	index, instack

begin
	# Allocate arrays used in depth first search

	call malloc (index, ADJ_SIZE(gp), TY_INT)
	call malloc (instack, ADJ_SIZE(gp), TY_INT)

	call amovkr (0.0, depth, ADJ_SIZE(gp))
	call amovki (NO, Memi[instack], ADJ_SIZE(gp))

	# If the starting node is zero, calculate the depth of every
	# node in the graph. Otherwise only calculate the depth of
	# nodes that are descendants of start

	if (start == 0) {
	    do root = 1, ADJ_SIZE(gp) {
   		if (depth[root] == 0.0)
	    	    call rootdepth (gp, root, Memi[index], Memi[instack], depth)
	    }
	} else {
	    call rootdepth (gp, start, Memi[index], Memi[instack], depth)
	}

	call mfree (instack, TY_INT)
	call mfree (index, TY_INT)

end

# ROOTDEPTH -- Calculate the depth of the descendants of root

procedure rootdepth (gp, root, index, instack, depth)

pointer	gp		#  i: Pointer to adjacency list structure
int	root		#  i: Root node for depth calculation
int	index[ARB]	# io: Index of current node in forward node list
int	instack[ARB]	# io: Flag indicating whether node is on stack
real	depth[ARB]	# io: Array containing node depths
#--
int	top, node1, node2
pointer	sp, errmsg, stack

string	incycle  "The following component is part of a cycle: %s %d"

begin
	# Allocate memory for error message

	call smark (sp)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	# Allocate and initialize stack used to hold nodes in a
	# depth first search

	call malloc (stack, ADJ_SIZE(gp), TY_INT)

	top = 0
	Memi[stack+top] = root

	# Push root on the stack

	depth[root] = 1.0
	index[root] = 1
	instack[root] = YES

	while (top >= 0) {

	    # Get the top node from the stack and its next descendant

	    node1 = Memi[stack+top]
	    node2 = ADJ_NODE(gp,node1,index[node1])

	    if (node2 != 0) {
		index[node1] = index[node1] + 1

		# Check for cycles in the graph

		if (instack[node2] == YES) {
		    call sprintf (Memc[errmsg], SZ_LINE, incycle)
		    call pargstr (ADJ_NAME(gp,node2))
		    call pargi (ADJ_NUMBER(gp,node2))
		    call error (1, Memc[errmsg])
		}

		# Push the node on the stack if depth is not yet maximized

		if (depth[node2] <= depth[node1]) {

		    top = top + 1
	 	    Memi[stack+top] = node2
		    depth[node2] = top + 1.0

		    index[node2] = 1
		    instack[node2] = YES
		}

	    # If no more descendants, pop top node from stack

	    } else {
		top = top - 1
		instack[node1] = NO
	    }
	}

	call mfree (stack, TY_INT)
	call sfree (sp)

end

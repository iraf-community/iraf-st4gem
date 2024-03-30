include	"../adjlist.h"

# NODEDIST -- Calculate the distance of nodes in a graph from the start
#
# The distance of a given node in a graph is the number of nodes that must be
# passed through from the starting node to reach the given node, inclusive
# of the given node.
#
# B.Simon	27-Jul-88	First Code

procedure nodedist (gp, start, dist)

pointer	gp		# i: Pointer to adjacency list structure
int	start		# i: Starting node for distance calculation
int	dist[ARB]	# o: Array containing node distances
#--
int	root
pointer	index

begin
	# Allocate array used in depth first search

	call malloc (index, ADJ_SIZE(gp), TY_INT)

	call amovki (0, dist, ADJ_SIZE(gp))

	# If the starting node is zero, calculate the distance of every
	# node in the graph. Otherwise only calculate the distance of
	# nodes that are descendants of start

	if (start == 0) {
	    do root = 1, ADJ_SIZE(gp) {
   		if (dist[root] == 0)
	    	    call rootdist (gp, root, Memi[index], dist)
	    }
	} else {
	    call rootdist (gp, start, Memi[index], dist)
	}

	call mfree (index, TY_INT)

end

# ROOTDIST -- Calculate the distance of the descendants of root

procedure rootdist (gp, root, index, dist)

pointer	gp		#  i: Pointer to adjacency list structure
int	root		#  i: Root node for dist calculation
int	index[ARB]	# io: Index of current node in forward node list
int	dist[ARB]	# io: Array containing node distances
#--
int	top, node1, node2
pointer	stack

begin
	# Allocate and initialize stack used to hold nodes in a
	# depth first search

	call malloc (stack, ADJ_SIZE(gp), TY_INT)

	top = 0
	Memi[stack+top] = root

	# Push root on the stack

	dist[root] = 1
	index[root] = 1

	while (top >= 0) {

	    # Get the top node from the stack and its next descendant

	    node1 = Memi[stack+top]
	    node2 = ADJ_NODE(gp,node1,index[node1])

	    if (node2 != 0) {
		index[node1] = index[node1] + 1

		# Push the node on the stack if it has not been visited

		if (dist[node2] == 0) {

		    top = top + 1
	 	    Memi[stack+top] = node2
		    dist[node2] = top + 1

		    index[node2] = 1
		}

	    # If no more descendants, pop top node from stack

	    } else {
		top = top - 1
	    }
	}

	call mfree (stack, TY_INT)

end

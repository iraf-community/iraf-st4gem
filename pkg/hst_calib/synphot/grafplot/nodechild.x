include	"../adjlist.h"

# NODECHILD -- Determine the children and spouses of nodes in a graph
#
# The child of a given node is any node which is a direct descendant of
# the given node. The spouse is any node which share one or more common
# children with the given node. Both the children and spouses are returned
# as arrays containing linked lists.
#
# B.Simon	09-Aug-88	First Code

procedure nodechild (gp, start, child, spouse)

pointer	gp		# i: Pointer to adjacency list structure
int	start		# i: Starting node for spousal search
int	child[ARB]	# o: Array containing linked list of children
int	spouse[ARB]	# o: Array containing linked list of spouses
#--
int	root
pointer	index, parent

begin
	# Allocate arrays used in depth first search

	call malloc (index, ADJ_SIZE(gp), TY_INT)
	call malloc (parent, ADJ_SIZE(gp), TY_INT)

	call amovki (0, child, ADJ_SIZE(gp))
	call amovki (0, spouse, ADJ_SIZE(gp))
	call amovki (0, Memi[parent], ADJ_SIZE(gp))

	# If the starting node is zero, calculate the spouse of every
	# node in the graph. Otherwise only calculate the spouse of
	# nodes that are descendants of start

	if (start == 0) {
	    do root = 1, ADJ_SIZE(gp) {
   		if (Memi[parent+root-1] == 0)
	    	    call rootchild (gp, root, Memi[index],
				    Memi[parent], child, spouse)
	    }
	} else {
	    call rootchild (gp, start, Memi[index],
			    Memi[parent], child, spouse)
	}

	call mfree (parent, TY_INT)
	call mfree (index, TY_INT)

end

# ROOTCHILD -- Calculate the children and spouses of the descendants of root

procedure rootchild (gp, root, index, parent, child, spouse)

pointer	gp		#  i: Pointer to adjacency list structure
int	root		#  i: Root node for spouse calculation
int	index[ARB]	# io: Index of current node in forward node list
int	parent[ARB]	# io: Array of node parents
int	child[ARB]	# io: Array containing linked list of children
int	spouse[ARB]	# io: Array of linked list node spouses
#--
int	top, node1, node2, node3
pointer	stack

begin

	# Allocate and initialize stack used to hold nodes in a
	# depth first search

	call malloc (stack, ADJ_SIZE(gp), TY_INT)

	top = 0
	Memi[stack+top] = root

	# Push root on the stack. Parent of root is special case

	index[root] = 1
	parent[root] = root

	while (top >= 0) {

	    # Get the top node from the stack and its next descendant

	    node1 = Memi[stack+top]
	    node2 = ADJ_NODE(gp,node1,index[node1])

	    if (node2 != 0) {
		index[node1] = index[node1] + 1

		# Push the node on the stack if not yet visited

		if (parent[node2] == 0) {

		    top = top + 1
	 	    Memi[stack+top] = node2

		    index[node2] = 1
		    parent[node2] = node1

		    # Insert node2 at the end of the child linked list
		    # of node1. The child linked list is linked to the
		    # of the spouse linked list and its starting link
		    # is flagged by being a negative number.

		    node3 = node1
		    while (spouse[node3] > 0)
			node3 = spouse[node3]

		    if (spouse[node3] == 0)
			spouse[node3] = - node2
		    else {
			node3 = -spouse[node3]
			while (child[node3] != 0)
			    node3 = child[node3]
			child[node3] = node2
		    }

		} else {

		    # This node has more than one parent
		    # Check to see if it is already in the
		    # linked list of spouses

		    node2 = parent[node2]
		    while (node1 != node2 && spouse[node2] > 0)
			node2 = spouse[node2]

		    # Link together the child and spouse linked lists
		    # of node1 and node2

		    if (node1 != node2) {

			node3 = node1
			while (spouse[node3] > 0)
			    node3 = spouse[node3]

			if (spouse[node3] == 0) {

			    spouse[node3] = spouse[node2]
			    spouse[node2] = node1

			} else {

			    node3 = - spouse[node3]
			    while (child[node3] != 0)
				node3 = child[node3]

			    child[node3] = - spouse[node2]
			    spouse[node2] = node1

			}
		    }
		}

	    # If no more descendants, pop top node from stack

	    } else {
		top = top - 1
	    }
	}

	call mfree (stack, TY_INT)

end

include	"../adjlist.h"

# NODEORDER -- Perform a topological sort on nodes in a graph
#
# The topological sort order of a given node in a graph is the reverse of
# the order in which that node is last visited by a depth first search.
#
# B.Simon	21-Jul-88	First Code

procedure nodeorder (gp, start, norder, order)

pointer	gp		# i: Pointer to adjacency list structure
int	start		# i: Starting node for order calculation
int	norder		# o: Length of output array
int	order[ARB]	# o: Array containing topological sort order
#--
int	root, count, front, back, temp
pointer	index, visit

begin
	# Allocate arrays used in depth first search

	call malloc (index, ADJ_SIZE(gp), TY_INT)
	call malloc (visit, ADJ_SIZE(gp), TY_INT)

	call amovki (0, order, ADJ_SIZE(gp))
	call amovki (NO, Memi[visit], ADJ_SIZE(gp))

	# If the starting node is zero, calculate the reverse sort 
	# order of every node in the graph. Otherwise only calculate 
	# the reverse sort order of nodes that are descendants of start.

	count = 0
	if (start == 0) {
	    do root = 1, ADJ_SIZE(gp) {
   		if (Memi[visit+root-1] == NO)
		    call rootorder (gp, root, count, Memi[index],
				    Memi[visit], order)
	    }
	} else {
	    call rootorder (gp, start, count, Memi[index], Memi[visit], order)
	}
        norder = count

	# Reverse the order array
	
	front = 1
	back = norder
	while (front < back) {
	    temp = order[front]
	    order[front] = order[back]
	    order[back] = temp
	    front = front + 1
	    back = back - 1
	}

	call mfree (index, TY_INT)
	call mfree (visit, TY_INT)

end

# ROOTORDER -- Calculate the reverse sort order of the descendants of root

procedure rootorder (gp, root, count, index, visit, order)

pointer	gp		#  i: Pointer to adjacency list structure
int	root		#  i: Root node for sort order calculation
int	count		#  i: Count of number of nodes already left
int	index[ARB]	# io: Index of current node in forward node list
int	visit[ARB]	# io: Flag indicating node already visited
int	order[ARB]	# io: Array containing reverse sort order
#--
int	top, node1, node2
pointer	stack

begin
	# Allocate and initialize stack used to hold nodes in a
	# depth first search
               
	call malloc (stack, ADJ_SIZE(gp), TY_INT)
	top = 0

	# Push root on the stack

	Memi[stack+top] = root
	index[root] = 1
	visit[root] = YES

	while (top >= 0) {

	    # Get the top node from the stack and its next descendant

	    node1 = Memi[stack+top]
	    node2 = ADJ_NODE(gp,node1,index[node1])

	    if (node2 != 0) {

		index[node1] = index[node1] + 1

		# Push the node on the stack if it has not been visited

		if (visit[node2] == NO) {

		    top = top + 1
	 	    Memi[stack+top] = node2

		    index[node2] = 1
		    visit[node2] = YES
		}

	    # If no more descendants, pop top node from stack
	    # and add to order array

	    } else {
		top = top - 1
		count = count + 1
		order[count] = node1
	    }
	}

	call mfree (stack, TY_INT)

end

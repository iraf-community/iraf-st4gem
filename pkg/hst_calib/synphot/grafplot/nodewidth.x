include	"../adjlist.h"

# NODEWIDTH -- Determine the width of nodes in a graph
#
# The width of a given node in a graph is its placement in the direction
# orthogonal to depth. Nodes are placed so that they are centered under
# their descendant nodes. If a node has no descendants, its width is one
# greater than the node with the greatest width already calculated at that
# depth.
#
# B.Simon	04-Aug-88	First Code

procedure nodewidth (gp, start, child, spouse, depth, width)

pointer	gp		# i: Pointer to adjacency list structure
int	start		# i: Starting node
int	child[ARB]	# i: Array containing linked list of children
int	spouse[ARB]	# i: Array containing linked list of spouses
real	depth[ARB]	# i: Array containing node depths
real	width[ARB]	# o: Array containing node widths
#--
int	root
pointer	index, visit
real	leftmost

begin
	# Allocate arrays used in depth first search

	call malloc (index, ADJ_SIZE(gp), TY_INT)
	call malloc (visit, ADJ_SIZE(gp), TY_INT)

	call amovkr (0.0, width, ADJ_SIZE(gp))
	call amovki (NO, Memi[visit], ADJ_SIZE(gp))

	# If the starting node is zero, calculate the width of every
	# node in the graph. Otherwise only calculate the width of
	# nodes that are descendants of start

	leftmost = 0.0
	if (start == 0) {
	    do root = 1, ADJ_SIZE(gp) {
   		if (Memi[visit+root-1] == 0)
	    	    call rootwidth (gp, root, child, spouse, depth,
				    Memi[index], Memi[visit], leftmost, width)
	    }
	} else {
	    call rootwidth (gp, start, child, spouse, depth,
			    Memi[index], Memi[visit], leftmost, width)
	}

	call mfree (visit, TY_INT)
	call mfree (index, TY_INT)

end

# ROOTWIDTH -- Calculate the width of the descendants of root

procedure rootwidth (gp, root, child, spouse, depth,
		     index, visit, leftmost, width)

pointer	gp		#  i: Pointer to adjacency list structure
int	root		#  i: Root node for spouse calculation
int	child[ARB]	#  i: Array containing linked list of children
int	spouse[ARB]	#  i: Linked list of node spouses
real	depth[ARB]	#  i: Array of node depths
int	index[ARB]	# io: Index of current node in forward node list
int	visit[ARB]	# io: Order in which nodes are first visited
real	leftmost	# io: Maximum width of any node
real	width[ARB]	#  o: Array of node widths
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

	index[root] = 1
	visit[root] = YES
             
	while (top >= 0) {

	    # Get the top node from the stack and its next descendant

	    node1 = Memi[stack+top]
	    node2 = ADJ_NODE(gp,node1,index[node1])

	    if (node2 != 0) {
		index[node1] = index[node1] + 1

		# Push the node on the stack if not yet visited

		if (visit[node2] == NO) {

		    top = top + 1
	 	    Memi[stack+top] = node2

		    index[node2] = 1
		    visit[node2] = YES

		}

	    } else {

		# Calculate width of all nodes in spouse linked list

		if (width[node1] == 0.0)
		    call calcwidth (node1, ADJ_SIZE(gp), child, spouse, 
				    depth, leftmost, width)

		# Pop node from top of stack

		top = top - 1
	    }
	}

	call mfree (stack, TY_INT)

end

# CALCWIDTH -- Calculate the width of all nodes in spouse linked list

procedure calcwidth (head, size, child, spouse, depth, leftmost, width)

int	head		#  i: Head node of the spouse linked list
int	size		#  i: Number of nodes in graph
int	child[ARB]	#  i: Array containing linked list of children
int	spouse[ARB]	#  i: Linked list of spouses
real	depth[ARB]	#  i: Array of node depths
real	leftmost	# io: Maximum width of any node
real	width[ARB]	# io: Array of node widths
#--
bool	gotminmax
int	node1, node2, nspouse, nmindepth, tail
real	mindepth, delta, minwidth, maxwidth, fstwidth, secwidth

begin
	# Compute the minimum depth of any node in the spouse chain

	mindepth = depth[head]
	for (node1 = spouse[head]; node1 > 0; node1 = spouse[node1])
	    mindepth = min (mindepth, depth[node1])

	# Count the number of nodes in the spouse chain

	nspouse = 0
	nmindepth = 0
	for (node1 = head; node1 > 0; node1 = spouse[node1]) {
	    nspouse = nspouse + 1
	    if (depth[node1] == mindepth)
		nmindepth = nmindepth + 1
	}

	# Find the minumum and maximum width among all children of
	# the head node

	tail = - node1
	gotminmax = false
	for (node1 = tail; node1 != 0; node1 = child[node1]) {
	    if (width[node1] != 0.0) {
		if (gotminmax) {
		    minwidth = min (minwidth, width[node1])
		    maxwidth = max (maxwidth, width[node1])
		} else {
		    gotminmax = true
		    minwidth = width[node1]
		    maxwidth = width[node1]
		}
	    }
	}

	if (! gotminmax) {
	    minwidth = leftmost + 1.0
	    maxwidth = leftmost + 1.0
	}

	# Center the spouse chain between the min and max width

	delta = (maxwidth - minwidth + 1.0) - real(nspouse)
	fstwidth = minwidth + max (0.0, 0.5 * delta)
	secwidth = max(fstwidth, maxwidth - real(nmindepth))

	# Set the widths of nodes in the spouse chain in order of increasing 
	# depth. This algorithm uses the comparison counting sort from 
	# Knuth's "Art of Computer Programming" Vol. 3, p. 76.

	for (node1 = head; node1 > 0; node1 = spouse[node1])
	    if (depth[node1] == mindepth)
		width[node1] = fstwidth
	    else
		width[node1] = secwidth


	for (node1 = spouse[head]; node1 > 0; node1 = spouse[node1])
	    for (node2 = head; node2 != node1; node2 = spouse[node2])
		if (depth[node1] < depth[node2])
		    width[node2] = width[node2] + 1.0
		else
		    width[node1] = width[node1] + 1.0

	# Update the leftmost width

	for (node1 = head; node1 > 0; node1 = spouse[node1])
	    leftmost = max (leftmost, width[node1])	

end

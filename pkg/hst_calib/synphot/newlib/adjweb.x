include	"../adjlist.h"

# ADJ_WEB -- Add the forward node list to the adjacency list
#
# Each node in the adjacency list  has associated with it a list of forward
# nodes which are its direct descendents. These nodes form the web of
# interconnections which are used to traverse the adjacency list.
#
# B. Simon	21-Jul-88	First Code

procedure adj_web (ngraph, innode, outnode, inidx, outidx, invert, gp)
		    
int	ngraph		#  i: Number of rows in graph
pointer	innode		#  i: Array of innode values
pointer	outnode		#  i: Array of outnode values
pointer	inidx		#  i: Index array sorted on innode and compname
pointer	outidx		#  i: Index array sorted on outnode
pointer	invert		#  i: Inverted index associating row with node
pointer	gp		# io: Pointer to a graph structure
#--
int	isort, osort, irow, orow, prow
int	oldnode, node, list

errchk	fib_malloc, fib_appendi

begin
	# Allocate memory for the forward node list

	do list = 1, ADJ_SIZE(gp)
	    call fib_malloc (ADJ_NODARY(gp,list), 1, TY_INT)

	# Fill in the for the forward node lists in the adjacency list

	irow = 1
	orow = 1
	prow = 1

	while (irow <= ngraph && orow <= ngraph) {
	    isort = Memi[inidx+irow-1]
	    osort = Memi[outidx+orow-1]

	    # Seach for two rows in the graph table where the input and
	    # output node numbers are the same

	    if (Memi[innode+isort-1] < Memi[outnode+osort-1]) {
		irow = irow + 1
		orow = prow

	    } else if (Memi[innode+isort-1] == Memi[outnode+osort-1]) {
		list = Memi[invert+osort-1]

		# Look for the end of the list of forward nodes

		oldnode = 0
		for (node = 1; ADJ_NODE(gp,list,node) > 0; node = node + 1)
		    oldnode = ADJ_NODE(gp,list,node)

		# Only add the node to the list if it is new

		if (oldnode != Memi[invert+isort-1])
		    call fib_appendi (Memi[invert+isort-1],
				      ADJ_NODARY(gp,list), node)

		if (orow < ngraph) {
		    orow = orow + 1
		} else {
		    irow = irow + 1
		    orow = prow
		}

	    } else {

		# Add a zero to the end of the list of forward nodes

		list = Memi[invert+osort-1]
		for (node = 1; ADJ_NODE(gp,list,node) > 0; node = node + 1)
		    ;
		call fib_appendi (0, ADJ_NODARY(gp,list), node)

		orow = orow + 1
		prow = orow
	    }
	}

	# Add a zero to the end of the remaining forward node lists

	while (orow <= ngraph) {
	    osort = Memi[outidx+orow-1]
	    list = Memi[invert+osort-1]

	    for (node = 1; ADJ_NODE(gp,list,node) > 0; node = node + 1)
		;

	    call fib_appendi (0, ADJ_NODARY(gp,list), node)
	    orow = orow + 1
	}
end

# ADJ_OPEN -- Open a graph table and read its contents into an adjacency list
#
# This procedure reads one or more graph tables into a data structure
# designed to simplify traversal of the graph.
#
# B.Simon	21-Jul-88	First Code

pointer procedure adj_open (graph)

char	graph[ARB]	# i: Name of the graph table
#--
int	ngraph, nkeyword, ncomp, inidx, outidx
pointer	innode, outnode, compname, keyword, gp

errchk	loadgraf, sortgraf, adj_build, mfree

begin
	# Load the graph table contents into arrays

	call loadgraf (graph, ngraph, nkeyword, ncomp, 
		       innode, outnode, compname, keyword)

	# Create sorted indices on input and ouput nodes

	call sortgraf (ngraph, ncomp, innode, outnode, compname,
		       inidx, outidx)

	# Build the adjacency list structure 

	call adj_build (ngraph, ncomp, compname, innode, inidx,
			outnode, outidx, gp)

	# Free dynamically allocated arrays

	call mfree (innode, TY_INT)
	call mfree (outnode, TY_INT)
	call mfree (compname, TY_CHAR)
	call mfree (keyword, TY_CHAR)
	call mfree (inidx, TY_INT)
	call mfree (outidx, TY_INT)

	return (gp)
end

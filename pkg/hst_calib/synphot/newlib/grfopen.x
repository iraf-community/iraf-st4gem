# GRF_OPEN -- Open a graph table and read its contents into a graph list
#
# This procedure reads one or more graph tables into an internal data structure
# whose structure parallels that of the graph table
#
# B.Simon	12-Aug-88	First Code

pointer procedure grf_open (graph)

char	graph[ARB]	# i: Name of the graph table
#--
int	ngraph, nkeyword, ncomp, inidx, outidx
pointer	innode, outnode, compname, keyword, gp

errchk	loadgraf, sortgraf, grf_build, mfree

begin
	# Load the graph table contents into arrays

	call loadgraf (graph, ngraph, nkeyword, ncomp, 
		       innode, outnode, compname, keyword)

	# Create sorted indices on input and ouput nodes

	call sortgraf (ngraph, ncomp, innode, outnode, compname,
		       inidx, outidx)

	# Build the graph list structure 

	call grf_build (ngraph, nkeyword, ncomp, innode, outnode, 
			compname, keyword, inidx, gp)

	# Free dynamically allocated arrays no longer needed

	call mfree (inidx, TY_INT)
	call mfree (outidx, TY_INT)

	return (gp)
end

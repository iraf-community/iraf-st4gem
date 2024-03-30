# SORTGRAF -- Create index arrays sorted on outnode, innode, and compname
#
# B. Simon	29-Jun-88	First Code

procedure sortgraf (ngraph, ncomp, innode, outnode, compname, inidx, outidx)

int	ngraph		# i: Number of rows in graph
int	ncomp		# o: Length of component name string
pointer	innode		# i: Array of innode values
pointer	outnode		# i: Array of outnode values
pointer	compname	# i: Array of component names
int	inidx		# o: Index array sorted on innode and compname
int	outidx		# o: Index array sorted on outnode
#--
int	igraph, mgraph

extern	int_compare(), str_compare()

errchk	malloc, mergesort

begin
	# Allocate dynamic memory to hold index arrays

	mgraph = 2 * ngraph
	call malloc (inidx, mgraph, TY_INT)
	call malloc (outidx, mgraph, TY_INT)

	# Sort the outidx array on outnode

	do igraph = 1, ngraph
	    Memi[outidx+igraph-1] = igraph

	call mergesort (Memi[outidx], outnode, 1, mgraph, ngraph, int_compare)
	call amovi (Memi[outidx], Memi[inidx], ngraph)

	# Sort the inidx array on innode, compname, and outnode

	call mergesort (Memi[inidx], compname,
			ncomp, mgraph, ngraph, str_compare)
	call mergesort (Memi[inidx], innode, 1, mgraph, ngraph, int_compare)

end

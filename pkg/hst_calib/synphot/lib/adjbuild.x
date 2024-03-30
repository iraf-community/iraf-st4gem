include	"../adjlist.h"

# ADJ_BUILD -- Build an adjacency list from the graph arrays
#
# B. Simon	21-Jul-88	First Code

procedure adj_build (ngraph, ncomp, compname, innode, inidx, 
		     outnode, outidx, gp)
		    
int	ngraph		# i: Number of rows in graph
int	ncomp		# i: Length of component name string
pointer	compname	# i: Array of component names
pointer	innode		# i: Array of innode values
pointer	inidx		# i: Index array sorted on innode and compname
pointer	outnode		# i: Array of outnode values
pointer	outidx		# i: Index array sorted on outnode
pointer	gp		# o: Pointer to the adjacency list structure
#--
int	igraph, isort, oldnode, maxlist
pointer	sp, oldcomp, invert

bool	strne()

errchk	salloc, calloc, malloc, mfree, adj_tags, adj_web

begin
	# Allocate memory to hold old component name

	call smark (sp)
	call salloc (oldcomp, ncomp, TY_CHAR)

	# Create the inverted index array

	call malloc (invert, ngraph, TY_INT)

	oldnode = 0
	maxlist = 0
	Memc[oldcomp] = EOS

	do igraph = 1, ngraph {

	    isort = Memi[inidx+igraph-1]

	    # Check to see if this is a new node

	    if (Memi[innode+(isort-1)] != oldnode ||
		strne (Memc[compname+(isort-1)*(ncomp+1)], Memc[oldcomp])) {

		oldnode = Memi[innode+(isort-1)]
		call strcpy (Memc[compname+(isort-1)*(ncomp+1)],
			     Memc[oldcomp], ncomp)
		maxlist = maxlist + 1
	    }

	    Memi[invert+(isort-1)] = maxlist
	}

	# Allocate memory to hold top level of the adjacency list

	gp = NULL
	call calloc (gp, LEN_ADJSTRUCT, TY_STRUCT)

	ADJ_SIZE(gp) = maxlist
	call malloc (ADJ_NODPTR(gp), maxlist, TY_INT)
	call malloc (ADJ_NUMPTR(gp), maxlist, TY_INT)
	call malloc (ADJ_NAMPTR(gp), maxlist, TY_INT)

	# Store the component names and node numbers in the adjacency list

	call adj_tags (ngraph, ncomp, compname, innode, inidx, invert, gp)

	# Fill in the list of forward nodes in the adjacency list

	call adj_web (ngraph, innode, outnode, inidx, outidx, invert, gp)

	# Free dynamic memory

	call mfree (invert, TY_INT)
	call sfree (sp)
end

include	"../adjlist.h"

# ADJ_TAGS -- Add the tag fields for the adjacency list structure
#
# This procedure creates the component name and node number fields associated 
# with each node in the instrument graph. These two fields are called tag
# fields because they are used to identify a node in the adjacency list
# structure.
#
# B. Simon	21-Jul-88	First Code

procedure adj_tags (ngraph, ncomp, compname, innode, inidx, invert, gp)
		    
int	ngraph		#  i: Number of rows in graph table
int	ncomp		#  i: Length of component name string
pointer	compname	#  i: Array of component names
pointer	innode		#  i: Array of innode values
pointer	inidx		#  i: Index array sorted on innode and compname
pointer	invert		#  i: Inverted index array associating row with node
pointer	gp		# io: Pointer to a graph structure
#--
int	oldnode, list, igraph, isort, reclen
pointer	record

int	strlen()

errchk	malloc

begin
	list = 0
	oldnode = 0
	do igraph = 1, ngraph {
	    isort = Memi[inidx+igraph-1]
	    
	    # Check to see if this is a new node

	    if (Memi[invert+(isort-1)] != oldnode) {
		
		# Update the old node number 

		oldnode = Memi[invert+(isort-1)]
		list = list + 1

		# Store the innode number in the adjacency list

		ADJ_NUMBER(gp,list) = Memi[innode+(isort-1)]

		# Store the component name in the adjacency list

		record = compname + (isort - 1) * (ncomp + 1)
		reclen = strlen (Memc[record])

		call malloc (ADJ_NAMARY(gp,list), reclen+1, TY_CHAR)
		call strcpy (Memc[record], ADJ_NAME(gp,list), reclen)
	    }
	}

end

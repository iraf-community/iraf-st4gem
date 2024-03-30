include	"../adjlist.h"

# WGRAFLIST -- Write the list of descendant graph nodes to STDOUT
#
# B.Simon	22-JUL-88	First Code

procedure wgraflist (gp, norder, order, dist)

pointer	gp		# i: Pointer to adjacency list
int	norder		# i: Number of nodes to write
int	order[ARB]	# i: Order in which nodes will be written
int	dist[ARB]	# i: Distance of nodes in the graph from start
#--
int	iorder, list, level

begin
	do iorder = 1, norder {
	    list = order[iorder]

	    # Indent the node according to its distance from start

      	    do level = 1, dist[list] - 1
		call putline (STDOUT, ". ")

	    # Write the component name

	    call putline (STDOUT, ADJ_NAME(gp,list))
	    call putci (STDOUT, '\n')
	}
end

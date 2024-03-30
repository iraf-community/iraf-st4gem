include	"../adjlist.h"

# GRAFLIST -- List all the descendants of a node in a graph table
#
# B.Simon	22-Jul-88	First Code

procedure graflist ()

pointer	grftable	# i: Graph table name
pointer compname	# i: Component name of starting node
#--
int	start, norder
pointer	sp, gp, order, dist

int	getnode()
pointer	adj_open()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (grftable, SZ_FNAME, TY_CHAR)
	call salloc (compname, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgstr ("grftable", Memc[grftable], SZ_FNAME)
	call clgstr ("compname", Memc[compname], SZ_FNAME)

	# Read the instrument graph into an adjacency list structure

	gp = adj_open (Memc[grftable])

	# Find the node whose descendants will be listed

	start = getnode (gp, Memc[compname])

	# Calculate the order in which the nodes will be listed
	# and the distance of each node to be listed

	call malloc (order, ADJ_SIZE(gp), TY_INT)
	call malloc (dist, ADJ_SIZE(gp), TY_INT)

	call nodeorder (gp, start, norder, Memi[order])
	call nodedist (gp, start, Memi[dist])

	# Write the list of descendant nodes to STDOUT

	call wgraflist (gp, norder, Memi[order], Memi[dist])

	# Free dynamic memory

	call mfree (dist, TY_INT)
	call mfree (order, TY_INT)
	call adj_close (gp)

	call sfree (sp)

end

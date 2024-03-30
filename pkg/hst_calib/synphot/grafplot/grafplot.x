include	"grafplot.h"
include	"../adjlist.h"

# GRAFPLOT -- Plot all the descendants of a node in a graph table
#
# B.Simon	22-Jul-88	First Code

procedure grafplot ()

pointer	grftable	# i: Graph table name
pointer compname	# i: Component name of starting node
pointer	title		# i: Plot title
pointer	device		# i: Graphics device name
#--
int	start
pointer	sp, gp, gd, width, depth, child, spouse
real	format[LENFMT]

int	getnode()
pointer	adj_open()
pointer	gopen()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (grftable, SZ_FNAME, TY_CHAR)
	call salloc (compname, SZ_FNAME, TY_CHAR)
	call salloc (title, SZ_FNAME, TY_CHAR)
	call salloc (device, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgstr ("grftable", Memc[grftable], SZ_FNAME)
	call clgstr ("compname", Memc[compname], SZ_FNAME)
	call clgstr ("title", Memc[title], SZ_FNAME)
	call clgstr ("device", Memc[device], SZ_FNAME)

	# Read the instrument graph into an adjacency list structure

	gp = adj_open (Memc[grftable])

	# Find the node whose descendants will be plotted

	start = getnode (gp, Memc[compname])

	# Calculate the location where the nodes will be plotted

	call malloc (width, ADJ_SIZE(gp), TY_REAL)
	call malloc (depth, ADJ_SIZE(gp), TY_REAL)
	call malloc (child, ADJ_SIZE(gp), TY_INT)
	call malloc (spouse, ADJ_SIZE(gp), TY_INT)

	call nodedepth (gp, start, Memr[depth])
	call nodechild (gp, start, Memi[child], Memi[spouse])
	call nodewidth (gp, start, Memi[child], Memi[spouse],
			Memr[depth], Memr[width])

	# Plot the descendant nodes

	gd = gopen (Memc[device], NEW_FILE, STDGRAPH)
	call gclear (gd)

	call plotscale (gd, gp, start, Memr[depth], Memr[width],
			Memc[title], format)
	call nodeplot (gd, format, gp, start, Memr[depth], Memr[width])

	call gclose (gd)

	# Free dynamic memory

	call mfree (spouse, TY_INT)
	call mfree (child, TY_INT)
	call mfree (depth, TY_REAL)
	call mfree (width, TY_REAL)

	call adj_close (gp)

	call sfree (sp)

end

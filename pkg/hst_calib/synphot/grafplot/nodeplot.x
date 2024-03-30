include	"../adjlist.h"
include	"grafplot.h"

# NODEPLOT -- Plot nodes in a graph
#
# B.Simon	05-Aug-88	First Code

procedure nodeplot (gd, format, gp, start, depth, width)

pointer	gd		# i: Graphics device structure
real	format[ARB]	# i: Array containing plot format information
pointer	gp		# i: Pointer to adjacency list structure
int	start		# i: Starting node for plot
real	depth[ARB]	# i: Array containing node depths
real	width[ARB]	# i: Array containing node widths
#--
int	root
pointer	index, plotted

begin
	# Allocate arrays used in depth first search

	call malloc (index, ADJ_SIZE(gp), TY_INT)
	call malloc (plotted, ADJ_SIZE(gp), TY_INT)

	call amovki (NO, Memi[plotted], ADJ_SIZE(gp))

	# If the starting node is zero, plot every node in the graph
	# Otherwise only plot the nodes that are descendants of start

	if (start == 0) {
	    do root = 1, ADJ_SIZE(gp) {
   		if (Memi[plotted+root-1] == NO)
	    	    call rootplot (gd, format, gp, root, depth, width,
				   Memi[index], Memi[plotted])
	    }
	} else {
	    call rootplot (gd, format, gp, start, depth, width,
			    Memi[index], Memi[plotted])
	}

	call mfree (plotted, TY_INT)
	call mfree (index, TY_INT)

end

# ROOTPLOT -- Plot the descendants of root

procedure rootplot (gd, format, gp, root, depth, width, index, plotted)

pointer	gd		#  i: Graphics device structure
real	format[ARB]	#  i: Array containing plot format information
pointer	gp		#  i: Pointer to adjacency list structure
int	root		#  i: Root node for plot
real	depth[ARB]	#  i: Array containing node depths
real	width[ARB]	#  i: Array containing node widths
int	index[ARB]	# io: Index of current node in forward node list
int	plotted[ARB]	# io: Flag indicating whether node has been plotted
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

	call nodebox (gd, format, gp, root, depth, width)
	plotted[root] = YES

	while (top >= 0) {

	    # Get the top node from the stack and its next descendant

	    node1 = Memi[stack+top]
	    node2 = ADJ_NODE(gp,node1,index[node1])

	    if (node2 != 0) {
		index[node1] = index[node1] + 1

	    	# Draw the component name inside a box

		if (plotted[node2] == NO) {

		    top = top + 1
	 	    Memi[stack+top] = node2

		    index[node2] = 1

		    call nodebox (gd, format, gp, node2, depth, width)
		    plotted[node2] = YES
		}

		# Draw a line between the parent and child node

		call nodeline (gd, format, node1, node2, depth, width)

	    # If no more descendants, pop top node from stack

	    } else {
		top = top - 1
	    }
	}

	call mfree (stack, TY_INT)

end

# NODEBOX -- Draw the node's component name inside a box

procedure nodebox (gd, format, gp, node, depth, width)

pointer	gd		# i: Graphics device structure
real	format[ARB]	# i: Array containing plot format information
pointer	gp		# i: Pointer to adjacency list structure
int	node		# i: Node to be plotted
real	depth[ARB]	# i: Array containing node depths
real	width[ARB]	# i: Array containing node widths
#--
pointer	sp, line1, line2, ic, jc, kc
real	xbox[5], ybox[5], xtext, ytext

begin
	# Allocate memory to hold component name

	call smark (sp)
	call salloc (line1, int(format[BDEPTH]), TY_CHAR)
	call salloc (line2, int(format[BDEPTH]), TY_CHAR)

	# Draw box around node

	xbox[1] = (format[BDEPTH] + format[SDEPTH]) * depth[node]
	xbox[2] = xbox[1] - format[BDEPTH]
	xbox[3] = xbox[2]
	xbox[4] = xbox[1]
	xbox[5] = xbox[1]

	ybox[1] = (format[BWIDTH] + format[SWIDTH]) * width[node]
	ybox[2] = ybox[1]
	ybox[3] = ybox[1] - format[BWIDTH]
	ybox[4] = ybox[3]
	ybox[5] = ybox[1]

	call gpline (gd, xbox, ybox, 5)

	# Break the component name in two at the last underscore
	# if it will not fit in the box

	kc = ADJ_NAMARY(gp,node)
	jc = kc
	for (ic = kc; Memc[ic] != EOS; ic = ic + 1)
	    if (Memc[ic] == '_')
		jc = ic + 1

	if ((ic - kc) < (int(format[BDEPTH]) - 2)) {
	    call strcpy (Memc[kc], Memc[line1], ic - kc)
	    Memc[line2] = EOS
	} else {
	    call strcpy (Memc[kc], Memc[line1], jc - kc)
	    call strcpy (Memc[jc], Memc[line2], ic - jc)
	}

	# Print the component name in the box

	xtext = (format[BDEPTH] + format[SDEPTH]) * depth[node] - 
		0.5 * format[BDEPTH]
	ytext = (format[BWIDTH] + format[SWIDTH]) * width[node] - 2.0

	call gtext (gd, xtext, ytext, Memc[line1], EOS)

	ytext = ytext - 2.0
	call gtext (gd, xtext, ytext, Memc[line2], EOS)

	# Free dynamic memory

	call sfree (sp)
end

# NODELINE -- Draw a line from node1 to node2

procedure nodeline (gd, format, node1, node2, depth, width)

pointer	gd		# i: Graphics device structure
real	format[ARB]	# i: Array containing plot format information
int	node1		# i: Node to be plotted
int	node2		# i: Node to be plotted
real	depth[ARB]	# i: Array containing node depths
real	width[ARB]	# i: Array containing node widths
#--
real	pos1[2], pos2[2]

begin
	pos1[1] = (format[BDEPTH] + format[SDEPTH]) * depth[node1]
	pos1[2] = (format[BWIDTH] + format[SWIDTH]) * width[node1] -
		  0.5 * format[BWIDTH]

	pos2[1] = (format[BDEPTH] + format[SDEPTH]) * depth[node2] -
		  format[BDEPTH]
	pos2[2] = (format[BWIDTH] + format[SWIDTH]) * width[node2] -
		  0.5 * format[BWIDTH]

	call gline (gd, pos1[1], pos1[2], pos2[1], pos2[2])

end

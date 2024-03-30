include	<gset.h>
include	"grafplot.h"
include	"../adjlist.h"

# PLOTSCALE -- Open the plotting device and set the scale
#
# B.Simon	05-Aug-88	First Code

procedure plotscale (gd, gp, start, depth, width, title, format)

pointer	gd		# i: Graphics device structure
pointer	gp		# i: Adjacency list structure
int	start		# i: Starting node of plot
real	depth[ARB]	# i: Array of node depths
real	width[ARB]	# i: Array of node widths
char	title[ARB]	# i: Plot title
real	format[ARB]	# o: Array containing plot format information
#--
int	before, after, list
pointer	ic
real	maxdepth, maxwidth, pgdepth, pgwidth, chdepth, chwidth
real	txsize, delta, x1, x2, y1, y2, xtext, ytext

real	ggetr()

begin
	# Compute the plot spacing and store in format array

	format[BDEPTH] = 4.0
	do list = 1, ADJ_SIZE(gp) {
	    if (width[list] > 0.0) {
		before = 0
		after = 0
		for (ic = ADJ_NAMARY(gp,list); Memc[ic] != EOS; ic = ic + 1) {
		    after = after + 1
		    if (Memc[ic] == '_') {
			before = before + after
			after = 0
		    }
		}
	        format[BDEPTH] = max (format[BDEPTH], 
				      real(before), real(after))
	    }
	}
	format[BDEPTH] = format[BDEPTH] + 2.0
	format[BWIDTH] = 6.0
	format[SDEPTH] = 0.5 * format[BDEPTH]
	format[SWIDTH] = 0.5 * format[BWIDTH]
	format[TWIDTH] = 3.0

	# Compute the page dimensions

	maxdepth = depth[1]
	maxwidth = width[1]
	do list = 2, ADJ_SIZE(gp) {
	    maxdepth = max (maxdepth, depth[list])
	    maxwidth = max (maxwidth, width[list])
	}

	pgdepth = (format[BDEPTH] + format[SDEPTH]) * maxdepth
	pgwidth = (format[BWIDTH] + format[SWIDTH]) * maxwidth + format[TWIDTH]

	# Compute the text character size

	chdepth = ggetr (gd, "ch")
	chwidth = ggetr (gd, "cw")

	txsize = min (1.0, 1.0 / (chdepth * pgdepth), 1.0 / (chwidth * pgwidth))
	txsize = int (txsize / TXINC) * TXINC
	if (txsize < TXMIN)
	    call error (1, "Too many components to fit on one page")

	# Set the text attributes

	if (txsize < 1.0)
	    call gseti (gd, G_TXQUALITY, GT_HIGH)

	call gseti (gd, G_TXHJUSTIFY, GT_CENTER)
	call gseti (gd, G_TXVJUSTIFY, GT_CENTER)
	call gsetr (gd, G_TXSIZE, txsize)

	# Set the viewport scale

	delta = (format[BWIDTH] + format[SWIDTH]) * width[start]
	x1 = 0.0
	x2 = 1.0 / (txsize * chdepth)
	y1 = delta - 1.0 / (2.0 * txsize * chwidth)
	y2 = delta + 1.0 / (2.0 * txsize * chwidth)

	call gswind (gd, x1, x2, y1, y2)

	# Write the title

	format[TWIDTH] = 3.0 / txsize

	if (title[1] != EOS) {
	    xtext = 0.5 * (x1 + x2)
	    ytext = 0.5 * format[TWIDTH] + y1

	    call gtext (gd, xtext, ytext, title, "q=n")
	}

end

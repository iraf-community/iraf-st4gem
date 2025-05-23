# Copyright restrictions apply - see st4gem$copyright.stsdas

include <error.h>
include <gset.h>
include	<mach.h>
include	"../atom.h"
include	"../neberr.h"
include	"ntplot.h" 

define	DEBUG	false

#--------------------------------------------------------------------20 Aug 97--
.help ntcontour.x Apr97 nebular/plot
.ih
NAME
.nf
  nt_contour - Interactive loop for diagnostic curves plots
 plot_update - Plot diagnostic curves over a specified range of ratios
.fi
.endhelp
#-------------------------------------------------------------------------------
#  NT_CONTOUR -- Interactive loop for diagnostic curves plots.

procedure nt_contour (pl, interactive)

#  Calling arguments:
pointer	pl			# I: plot state descriptor
bool	interactive		# I: interactive session?

#  Declarations:
int	clgcur()		# CL graphics cursor
char	command[SZ_LINE]	# input cursor command
pointer	gp			# graphics descriptor
int	key			# input cursor key 
pointer	o			# parameter object
pointer	pr_alloc()		# allocate a parameter object
int	wcs			# index of WCS coordinate system
real	x_pos, y_pos		# cursor X, Y world coordinates

errchk	plot_setup

begin
	call plot_setup (pl)
	if (DEBUG)
	    call plt_debug (pl)
	call plot_update (pl)

	if (!interactive) 
	    return

	gp = PL_GP(pl)

	while (clgcur ("cur", x_pos, y_pos, wcs, key, command, SZ_LINE) != EOF) {
	    switch (key) {

	    # Report cursor position in world coordinates
	    case 'c':
		call printf ("X: %6.3g; Y: %6.3g\n")
		call pargr (x_pos)
		call pargr (y_pos)

	    # Quit task
	    case 'q', 'I':
		break

	    # Refresh plot
	    case 'r':
		call gframe (gp)
		call plot_setup (pl)
		iferr (call plot_update (pl)) {
		    call erract (EA_WARN)
		}

	    # Print help for keystroke/colon commands
	    case '?':
		call gdeactivate (gp, 0)
		call gpagefile (gp, KEYHELP, PROMPT)
		call greactivate (gp, 0)

	    # Colon command
	    case ':':
		o = pr_alloc (false)
		call strlwr (command)
		call pr_sset (o, command)
		iferr (call nt_colon (pl, o)) {
		    call erract (EA_WARN)
		}
		call pr_free (o)

	    # Unrecognized cursor command -- ring bell.
	    default:
		call printf ("Invalid cursor command\7\n")
	    }
	}
end


#-------------------------------------------------------------------------------
#  PLOT_UPDATE -- Plot diagnostic curves over a specified range of ratios

procedure plot_update (pl)

#  Calling arguments:
pointer	pl			# I: plot state descriptor

#  Declarations:
int	curve_calc()		# calculate curve for the diagnostic ratio
real	dratio			# delta diagnostic line ratio
pointer	gp			# graphics descriptor
int	i			# loop index
int	n_curve			# no. secondary contours to plot
int	n_pts			# no. points calculated per curve
int	n_plot			# no. points actually plotted @curve
int	nscan()			# no. tokens successfully converted

errchk	curve_calc, curve_plot

begin
	gp    = PL_GP(pl)
	n_pts = PL_NPTS(pl)
	call gseti (gp, G_PLTYPE, GL_SOLID)

	# Update labels & title
	if (!PL_APPEND(pl)) {
	    call label_update (pl)
	    call title_update (pl)
	    call glabax (gp, TITLE(pl), XLABEL(pl), YLABEL(pl))
	}

	# Plot reference curve. 
#	call printf ("Diag. expression: %s")
#	    call pargstr (DIAG_EXPR(pl))
	n_plot = curve_calc (pl, PL_RATIO(pl), DIAG_EXPR(pl))

	# Error out if no valid points in reference curve were plotted.  
	if (n_plot < 1)
	    call error (INVLD_PLOT, 
		"No valid points to plot: check input ratio or N_e/T_e limits")
	else
	    call curve_plot (pl, GL_DOTTED, PL_REF_COLOR(pl))

	# Plot the surrounding secondary curves. 
	call sscan (PL_LIST(pl))
	call gargr (dratio)
	n_curve = 0
	i = nscan ()
	while (i > n_curve) {
	    n_plot = curve_calc (pl, dratio, DIAG_EXPR(pl))
	    call curve_plot (pl, GL_SOLID, PL_DELTA_COLOR(pl))
	    n_curve = n_curve + 1
	    call gargr (dratio)
	    i = nscan ()
	} 

	# Clean up.
	call gflush (gp)
end



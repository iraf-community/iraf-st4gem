# Copyright restrictions apply - see st4gem$copyright.stsdas

include <gset.h>
include	<tbset.h>
include	"../atom.h"
include	"../fivel.h"
include	"../neberr.h"
include	"ntplot.h" 
include	"../zones.h"

#---------------------------------------------------------------------7 Sep 97--
.help diagplot.x Feb97 nebular/plot
.ih
NAME
.nf
   diagplot - Loop over input diagnostic ratios & plot
.fi
.endhelp
#-------------------------------------------------------------------------------
#  DIAGPLOT --	Loop over input diagnostic ratios & plot

procedure diagplot (ilr, otp, name)

#  Calling arguments:
pointer	ilr			# input line ratio descriptor
pointer	otp			# output table descriptor
char	name[ARB]		# specific object name

#  Local variables:
int	atom			# atomic number
int	ion			# spectrum (ion) number
int	clgeti()		# fetch TY_INT value from the CL
char	colname[SZ_COLNAME]	# output table column name
pointer	colptr[5]		# reference table column pointers
int	curve_calc()		# calculate curve for the diagnostic ratio
char	device[SZ_FNAME]	# graphics device
real	diag_ratio		# diagnostic ratio
#int	diag_type		# type of diagnostic 
char	expr[SZ_LINE]		# diagnostic expression
real	get_dratio()		# fetch diagnostic ratio from ref table
pointer gopen()                 # Returns pointer to graphics structure
pointer	gp			# graphics descriptor
int	i			# generic
int	line_type		# line plot type
int	n_pts			# no. points calculated per curve
int	n_plot			# no. points actually plotted @ curve
int	n_rows			# number of rows in reference table 
pointer	pl			# local plot descriptor
pointer	plt_alloc()		# allocate plot descriptor
char    title[SZ_LINE]		# Title for plot
pointer	tp			# reference curve table descriptor
int	zone			# ionization zone

define	MIN_PTS		10	# Minimum useful curve resolution

errchk	curve_calc, curve_plot, curve_to_table, gopen, open_curve, plot_setup 

begin
	n_pts = clgeti ("resolution")
	if (n_pts < MIN_PTS) 
	    call error (BAD_PLOT_RESOL, "Insufficient resolution")

	# Initialize graphics attributes
	pl = plt_alloc (n_pts)
	PL_TYPE(pl) = TE_VS_NE
	call clgstr ("device", device, SZ_FNAME)
	gp = gopen (device, NEW_FILE, STDGRAPH)
	PL_GP(pl) = gp

        call gseti (gp, G_NTITLELINES, 1)
        call gseti (gp, G_PLTYPE, GL_SOLID)

	# Fetch desired plot attributes from the CL and compute 
	# reference density & temperature arrays. 
	call get_basic_pltpars (pl, false)
	call plot_setup (pl)

	# Set label attributes. 
	call strcpy ("Nebular Diagnostic Diagram\n for: ", title,  SZ_LINE)
	call strcat (name, title, SZ_LINE)
	call strcpy ("Electron Density", XLABEL(pl), SZ_LINE)
	call strcpy ("Electron Temperature", YLABEL(pl), SZ_LINE)
	call glabax (gp, title, XLABEL(pl), YLABEL(pl))

	# Write reference density, temperature to output table. 
	call curve_to_table (otp, "Ref_Ne", DENSITY, NE_REF(pl), n_pts)
	call curve_to_table (otp, "Ref_Te", TEMPERATURE, TE_REF(pl), n_pts)

	# Plot density curves.
	call strcpy (EOS, expr, SZ_LINE)
	call open_curve (tp, colptr, n_rows)

	do i = 1, n_rows {
	    call get_diag_line (tp, colptr, i, atom, ion, zone, 
		PL_DIAG_TYPE(pl), colname)

	    diag_ratio = get_dratio (ilr, atom, ion, PL_DIAG_TYPE(pl))
	    if ( IS_INDEFR (diag_ratio) )
		next

	    if (zone == LOW)
	    	line_type = GL_SOLID

	    else if (zone == MED)
	    	line_type = GL_DASHED

	    else if (zone == HI)
	    	line_type = GL_DOTTED

	    call at_alloc (atom, ion, PL_AT(pl))
	    n_plot = curve_calc (pl, diag_ratio, expr)
	    if (n_plot < 1)
		next

	    call curve_plot (pl, line_type, 1)

	    if (PL_DIAG_TYPE(pl) == DENSITY)
	    	call curve_to_table (otp, colname, PL_DIAG_TYPE(pl), PL_NE(pl), 
					n_pts)

	    else if (PL_DIAG_TYPE(pl) == TEMPERATURE)
	    	call curve_to_table (otp, colname, PL_DIAG_TYPE(pl), PL_TE(pl), 
					n_pts)
	}

	# Clean up.
	call tbtclo (tp)
	call gflush (gp)
	call plt_free (pl)
end



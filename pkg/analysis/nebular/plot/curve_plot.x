# Copyright restrictions apply - see st4gem$copyright.stsdas

include	"../atom.h" 
include	"ntplot.h" 

define	VERB	false
define	TOLER	2.e-4

#---------------------------------------------------------------------6 Jul 97--
.help curve_plot.x Mar97 nebular/plot
.ih
NAME
.nf
    curve_calc - Compute density, temperature, or intensity array
    curve_plot - Plot diagnostic curve
.fi
.endhelp
#-------------------------------------------------------------------------------
#  CURVE_CALC -	Compute dens/temper array & plot buffer.  Returns number of 
#		points successfully computed.  

int procedure curve_calc (pl, diag_ratio, diag_expr)

#  Calling arguments:
pointer	pl			# I: plot descriptor
real	diag_ratio		# I: diagnostic line ratio
char	diag_expr[ARB]		# I: expression for diagnostic ratio
int	n_curve			# O: no. non-INDEF points in computed curve

#  Declarations:
pointer	at			# atomic data object
int	calc_density()		# routine to calculate density array
int	calc_intensity()	# routine to calculate intensity array
int	calc_temperature()	# routine to calculate temperature array
int	npts			# size of arrays
pointer	ref			# reference N_e/T_e array
pointer	sp			# top of stack memory
bool	streq()			# are two strings equal?
char	sx[SZ_LINE]		# generic

define	Ref	Memr[ref]

errchk	calc_density, calc_intensity, calc_temperature

begin
	n_curve = 0
	if ( IS_INDEFR (diag_ratio) )
	    return (n_curve)

	npts = PL_NPTS(pl)

	at = PL_AT(pl)
	if (streq (EOS, diag_expr) )
	    call set_diag_expr (at, PL_DIAG_TYPE(pl), sx, SZ_LINE)

	else
	    call strcpy (diag_expr, sx, SZ_LINE)

	if (PL_DIAG_TYPE(pl) == DENSITY)
	    n_curve = calc_density (at, TE_REF(pl), PL_NE(pl), npts, 
					diag_ratio, sx, TOLER, VERB)

	else if (PL_DIAG_TYPE(pl) == TEMPERATURE)
	    n_curve = calc_temperature (at, NE_REF(pl), PL_TE(pl), npts, 
					diag_ratio, sx, TOLER, VERB)

	else if (PL_DIAG_TYPE(pl) == INTENSITY) {
	    call smark (sp)
	    call salloc (ref, npts, TY_REAL)
	    call amovkr (diag_ratio, Ref, npts)

	    if (PL_TYPE(pl) == INTENS_VS_NE) 
	    	n_curve = calc_intensity (at, NE_REF(pl), Ref, PL_INT(pl), 
					npts, sx)

	    else if (PL_TYPE(pl) == INTENS_VS_TE) 
	    	n_curve = calc_intensity (at, Ref, TE_REF(pl), PL_INT(pl), 
					npts, sx)

	    call sfree (sp)
	}

	return (n_curve)
end


#-------------------------------------------------------------------------------
#  CURVE_PLOT -	Determines the applicable arrays & plots. 

procedure curve_plot (pl, line_type, line_color)

#  Calling arguments:
pointer	pl			# I: plot descriptor
int	line_type		# I: 
int	line_color		# I: 

# Declarations:
pointer	gp			# graphics descriptor
pointer	x_arr, y_arr		# X,Y curve arrays

# Memory management
define	X_arr	Memr[x_arr]
define	Y_arr	Memr[y_arr]

include <gset.h>

begin
	gp = PL_GP(pl)

	if (PL_DIAG_TYPE(pl) == DENSITY) {
	    x_arr = NE_PTR(pl)
	    y_arr = TE_REF_PTR(pl)

	} else if (PL_DIAG_TYPE(pl) == TEMPERATURE) {
	    x_arr = NE_REF_PTR(pl)
	    y_arr = TE_PTR(pl)

	} else if (PL_DIAG_TYPE(pl) == INTENSITY) {
	    x_arr = NE_REF_PTR(pl)
	    y_arr = INT_PTR(pl)

	    if (PL_TYPE(pl) == INTENS_VS_TE) 
	    	x_arr = TE_REF_PTR(pl)
	}

	call gseti (gp, G_PLTYPE, line_type)
	call gseti (gp, G_PLCOLOR, line_color)
	call gpline (gp, X_arr, Y_arr, PL_NPTS(pl))
	call gflush (gp)

	# Reset to default graphics types
	call gseti (gp, G_PLTYPE, GL_SOLID)
	call gseti (gp, G_PLCOLOR, 1)
end



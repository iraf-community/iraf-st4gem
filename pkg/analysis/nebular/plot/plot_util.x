# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<evvexpr.h>
include	"../at.h" 
include	"../neberr.h" 
include	"ntplot.h" 

#---------------------------------------------------------------------9 Jul 97--
.help plot_util.x Mar97 nebular/plot
.ih
NAME
.nf
    plot_setup - Set plot attributes and initialize reference arrays
      p_update - Update plot reference array
   p_eval_expr - Evaluate expression for line flux ratio
  ratio_update - Update the diagnostic line ratio
   diag_update - Update the diagnostic expression
  label_update - Update the plot labels for ntplot display
  title_update - Update the plot title for ntplot display
      mk_title - Construct default plot title for ntplot display
.fi
.endhelp
#-------------------------------------------------------------------------------
#  PLOT_SETUP -- Set plot attributes and initialize reference arrays. 

procedure plot_setup (pl)

include <gset.h>

#  Arguments:
pointer	pl			# I: plot descriptor

#  Declarations:
pointer	gp			# graphics descriptor
int	n_pts			# size of arrays
bool	plot_log_x		# plot log X axis?
bool	plot_log_y		# plot log Y axis?

begin
	n_pts = PL_NPTS(pl)

	# Compute reference densities, temperatures & intensities
	switch (PL_TYPE(pl)) {
	case TE_VS_NE:
	    plot_log_x = PLOT_LOG_NE(pl)
	    plot_log_y = PLOT_LOG_TE(pl)

	    call p_update (NE_REF(pl), n_pts, plot_log_x, LO_NE(pl), HI_NE(pl))
	    call p_update (TE_REF(pl), n_pts, plot_log_y, LO_TE(pl), HI_TE(pl))

	case INTENS_VS_TE:
	    plot_log_x = PLOT_LOG_TE(pl)
	    plot_log_y = PLOT_LOG_INTENS(pl)

	    call p_update (TE_REF(pl), n_pts, plot_log_x, LO_TE(pl), HI_TE(pl))
	    call p_update (INT_REF(pl), n_pts, plot_log_y, LO_INTENS(pl), 
			HI_INTENS(pl))

	case INTENS_VS_NE:
	    plot_log_x = PLOT_LOG_NE(pl)
	    plot_log_y = PLOT_LOG_INTENS(pl)

	    call p_update (NE_REF(pl), n_pts, plot_log_x, LO_NE(pl), HI_NE(pl))
	    call p_update (INT_REF(pl), n_pts, plot_log_y, LO_INTENS(pl), 
			HI_INTENS(pl))
	}

	gp = PL_GP(pl)
	if (plot_log_x) 
	    call gseti (gp, G_XTRAN, GW_LOG)
	else 
            call gseti (gp, G_XTRAN, GW_LINEAR)

	if (plot_log_y) 
	    call gseti (gp, G_YTRAN, GW_LOG)
	else 
            call gseti (gp, G_YTRAN, GW_LINEAR)

	if (!PL_APPEND(pl)) {
	    switch (PL_TYPE(pl)) {
	    case TE_VS_NE:
	    	call gswind (gp, LO_NE(pl), HI_NE(pl), LO_TE(pl), HI_TE(pl))

	    case INTENS_VS_NE:
	    	call gswind (gp, LO_NE(pl), HI_NE(pl), LO_INTENS(pl), 
				HI_INTENS(pl))

	    case INTENS_VS_TE:
	    	call gswind (gp, LO_TE(pl), HI_TE(pl), LO_INTENS(pl), 
				HI_INTENS(pl))
	    }
	}
end


#-------------------------------------------------------------------------------
#  P_UPDATE -- Update plot reference array. 

procedure p_update (array, n_pts, plot_log, x1, x2)

include <gset.h>

#  Arguments:
real	array[ARB]		# reference array
int	n_pts			# size of arrays
bool	plot_log		# plot axis in log?
real	x1, x2			# upper/lower limits for array

#  Declarations:
int	i			# generic
real	incr			# increment for array values
real	log_x1, log_x2		# log10 of x1, x2

begin
	# Even steps in log
	if (plot_log) {
	    log_x1 = log10 (x1)
	    log_x2 = log10 (x2)
	    incr = (log_x2 - log_x1) / (n_pts - 1.)
	    do i = 1, n_pts
	    	array[i] = 10 ** (log_x1 + (i - 1) * incr)

	} else {
	    incr = (x2 - x1) / (n_pts - 1)
	    do i = 1, n_pts
	    	array[i] = x1 + (i - 1) * incr
	}
end


#-------------------------------------------------------------------------------
#  P_EVAL_EXPR - Evaluate expression for line flux ratio. 

real procedure p_eval_expr (expression)

include	<evvexpr.h>

#  Arguments:
char	expression[ARB]		# I: expression to evaluate
real	ratio			# O: diagnostic line ratio

#  Declarations:
pointer	ep			# pointer to expression structure
pointer	evvexpr()		# evaluate an algebraic expression

errchk	evvexpr

begin
	ratio = INDEFR
	ep = evvexpr (expression, NULL, NULL, NULL, NULL, EV_RNGCHK)

	switch (O_TYPE(ep)) {
	case TY_REAL:
	    ratio = O_VALR(ep)

	default:
	    call error (BAD_RATIO, 
		"Expression must involve floating-point numbers")
	}

	call mfree (ep, TY_STRUCT)
	return (ratio)
end


#-------------------------------------------------------------------------------
#  RATIO_UPDATE - Update the diagnostic line ratio.

procedure ratio_update (pl, expstr)

#  Arguments:
pointer	pl			# I: plot state descriptor
char	expstr[ARB]		# I: expression to evaluate

#  Declarations:
real	p_eval_expr ()		# evaluate an expression

begin
	iferr ( PL_RATIO(pl) = p_eval_expr (expstr) ) 
	    call error (BAD_RATIO, "Invalid expression for reference curve")

	if (IS_INDEFR(PL_RATIO(pl))) 
	    call error (BAD_RATIO, "Must provide diagnostic line ratio")
end


#-------------------------------------------------------------------------------
#  DIAG_UPDATE - Update the diagnostic expression.

procedure diag_update (pl, expression, max_char)

#  Arguments:
pointer	pl			# I: plot state descriptor
char	expression[max_char]	# I/O: expression to update
int	max_char		# I: max chars in expression

#  Declarations:
int	dtype			# type of diagnostic
bool	streq()			# are two strings equal?

begin
	call strupr (expression)
	if ( streq (expression, "DEFAULT") || streq (expression, EOS) ) {
	    if (PL_TYPE(pl) == INTENS_VS_NE)
		dtype = DENSITY

	    else if (PL_TYPE(pl) == INTENS_VS_TE)
		dtype = TEMPERATURE

	    else 
		dtype = PL_DIAG_TYPE(pl)

	    call set_diag_expr (PL_AT(pl), dtype, DIAG_EXPR(pl), max_char)

	} else 
	    call strcpy (expression, DIAG_EXPR(pl), max_char)
end


#-------------------------------------------------------------------------------
#  LABEL_UPDATE - Update the plot labels for ntplot display.

procedure label_update (pl)

#  Argument:
pointer	pl			# I: plot state descriptor

begin
	if (PL_TYPE(pl) == TE_VS_NE) {
	    call strcpy ("Electron Density", XLABEL(pl), SZ_LINE)
	    call strcpy ("Electron Temperature", YLABEL(pl), SZ_LINE)

	} else if (PL_TYPE(pl) == INTENS_VS_NE) {
	    call strcpy ("Electron Density", XLABEL(pl), SZ_LINE)
	    call strcpy ("Intensity Ratio", YLABEL(pl), SZ_LINE)

	} else if (PL_TYPE(pl) == INTENS_VS_TE) {
	    call strcpy ("Electron Temperature", XLABEL(pl), SZ_LINE)
	    call strcpy ("Intensity Ratio", YLABEL(pl), SZ_LINE)
	}
end


#-------------------------------------------------------------------------------
#  TITLE_UPDATE - Update the plot title for ntplot display.

procedure title_update (pl)

#  Arguments:
pointer	pl			# I: plot state descriptor

# Declarations:
pointer	at			# atomic data object
int	get_atom_str()		# Get the atom/ion string
int	get_ion_str()		# Get the ion diagnostic string
char	ion_str[SZ_LINE]	# specific transition description
int	nchars			# no. characters in returned string
bool	streq()			# are two strings equal?
char	sx[SZ_LINE]		# generic

begin
	if ( streq (USER_TRANSITION(pl), "DEFAULT") || 
	     streq (USER_TRANSITION(pl), EOS) ) {
	    if (PL_TYPE(pl) == INTENS_VS_NE)
	    	nchars = get_ion_str (PL_AT(pl), DENSITY, ion_str, SZ_LINE)

	    else if (PL_TYPE(pl) == INTENS_VS_TE)
	    	nchars = get_ion_str (PL_AT(pl), TEMPERATURE, ion_str, SZ_LINE)

	    else 
	    	nchars = get_ion_str (PL_AT(pl), PL_DIAG_TYPE(pl), ion_str, 
					SZ_LINE)
	} else {
	    at = PL_AT(pl)
	    nchars = get_atom_str (AT_ATOM(at), AT_ION(at), sx, SZ_LINE)
	    call sprintf (ion_str, SZ_LINE, "%s: %s")
		call pargstr (sx)
		call pargstr (DIAG_EXPR(pl))
	}

	call mk_title (PL_RATIO(pl), PL_LIST(pl), ion_str, USER_TITLE(pl), 
			TITLE(pl))
end


#-------------------------------------------------------------------------------
#  MK_TITLE - Construct default plot title for ntplot display.

procedure mk_title (ratio, list, ion_str, title2, title)

#  Calling arguments:
real	ratio			# I: diagnostic line ratio
char	list[ARB]		# I: list of secondary contours
char	ion_str[SZ_LINE]	# I: specific transition
char    title2[SZ_LINE]		# I: second line (user title) for plot
char    title[3*SZ_LINE]	# O: Title for plot

#  Declarations:
bool	streq()			# are two strings equal?
char    title1[SZ_LINE]		# first line of title for plot
char    title3[SZ_LINE]		# third line of title for plot

include	<mach.h>

begin
	# Construct title string.
	# First line: 
	call strcpy ("Diagnostic Ratio Contour Plot", title1, SZ_LINE)

	# Third line: 
	call sprintf (title3, SZ_LINE, "Primary = %g; Secondary = %s")
	    call pargr (ratio)
	if (streq (list, EOS) ) 
	    call pargstr ("[none]")
	else
	    call pargstr (list)

	# String 'em together.
	call sprintf (title, SZ_LINE, "%s\n %s\n %s")
	call pargstr (title1)

	# Construct default diagnostic expression...
	if (streq (title2, "default") || streq (title2, EOS) ) 
	    call pargstr (ion_str)

	# ...or user expression.
	else 
	    call pargstr (title2)

	call pargstr (title3)
end



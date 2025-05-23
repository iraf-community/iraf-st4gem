# Copyright restrictions apply - see st4gem$copyright.stsdas

include	"ntplot.h" 

#---------------------------------------------------------------------9 Jul 97--
.help plot_struct.x Mar97 nebular/plot
.ih
NAME
.nf
   plt_alloc - Allocate plot data structure
 plt_realloc - Re-allocate reference array sizes in plot data structure
    plt_free - De-allocate plot data structure
.fi
.ih
DESCRIPTION
.endhelp
#-------------------------------------------------------------------------------
#  PLT_ALLOC - Allocate plot data structure. 

pointer procedure plt_alloc (n_pts)

#  Arguments:
int	n_pts			# I: no. points calculated per curve
pointer	pl			# O: plot descriptor

begin
	call calloc (pl, LEN_PLT, TY_STRUCT)

	PL_TYPE(pl) = INDEFI
	PL_APPEND(pl) = false
	PL_GP(pl) = NULL
	if ( IS_INDEFI(n_pts) || n_pts < MIN_PTS)
	    PL_NPTS(pl) = MIN_PTS
	else
	    PL_NPTS(pl) = n_pts

	call calloc (NE_PTR(pl),      PL_NPTS(pl), TY_REAL)
	call calloc (TE_PTR(pl),      PL_NPTS(pl), TY_REAL)
	call calloc (INT_PTR(pl),     PL_NPTS(pl), TY_REAL)
	call calloc (NE_REF_PTR(pl),  PL_NPTS(pl), TY_REAL)
	call calloc (TE_REF_PTR(pl),  PL_NPTS(pl), TY_REAL)
	call calloc (INT_REF_PTR(pl), PL_NPTS(pl), TY_REAL)
	call calloc (PL_LIST_PTR(pl),    SZ_LINE, TY_CHAR)
	call calloc (X_LABEL_PTR(pl),    SZ_LINE, TY_CHAR)
	call calloc (Y_LABEL_PTR(pl),    SZ_LINE, TY_CHAR)
	call calloc (TITLE_PTR(pl),    3*SZ_LINE, TY_CHAR)
	call calloc (USER_TITLE_PTR(pl), SZ_LINE, TY_CHAR)
	call calloc (USER_TRANS_PTR(pl), SZ_LINE, TY_CHAR)
	call calloc (DIAG_EXPR_PTR(pl),  SZ_LINE, TY_CHAR)
	call calloc (PL_TABLE_PTR(pl),   SZ_FNAME, TY_CHAR)

	return (pl)
end


#-------------------------------------------------------------------------------
#  PLT_REALLOC - Re-allocate reference array sizes in plot data structure. 

procedure plt_realloc (pl, n_pts)

#  Arguments:
pointer	pl			# I: plot descriptor
int	n_pts			# I: no. points calculated per curve

begin
	if ( IS_INDEFI(n_pts) || n_pts < MIN_PTS)
	    PL_NPTS(pl) = MIN_PTS
	else
	    PL_NPTS(pl) = n_pts

	call realloc (NE_PTR(pl),      PL_NPTS(pl), TY_REAL)
	call realloc (TE_PTR(pl),      PL_NPTS(pl), TY_REAL)
	call realloc (INT_PTR(pl),     PL_NPTS(pl), TY_REAL)
	call realloc (NE_REF_PTR(pl),  PL_NPTS(pl), TY_REAL)
	call realloc (TE_REF_PTR(pl),  PL_NPTS(pl), TY_REAL)
	call realloc (INT_REF_PTR(pl), PL_NPTS(pl), TY_REAL)
end


#-------------------------------------------------------------------------------
#  PLT_FREE - De-allocate plot data structure. 

procedure plt_free (pl)

#  Arguments:
pointer	pl			# I: plot descriptor

begin
	call mfree (NE_PTR(pl), TY_REAL)
	call mfree (TE_PTR(pl), TY_REAL)
	call mfree (INT_PTR(pl), TY_REAL)
	call mfree (NE_REF_PTR(pl), TY_REAL)
	call mfree (TE_REF_PTR(pl), TY_REAL)
	call mfree (INT_REF_PTR(pl), TY_REAL)

	call mfree (PL_LIST_PTR(pl),    TY_CHAR)
	call mfree (X_LABEL_PTR(pl),    TY_CHAR)
	call mfree (Y_LABEL_PTR(pl),    TY_CHAR)
	call mfree (TITLE_PTR(pl),      TY_CHAR)
	call mfree (USER_TITLE_PTR(pl), TY_CHAR)
	call mfree (USER_TRANS_PTR(pl), TY_CHAR)
	call mfree (DIAG_EXPR_PTR(pl),  TY_CHAR)
	call mfree (PL_TABLE_PTR(pl),   TY_CHAR)

	if (PL_GP(pl) != NULL)
	    call gclose (PL_GP(pl))

	if (PL_AT(pl) != NULL)
	    call at_free (PL_AT(pl))

	call mfree (pl, TY_STRUCT)
end


#-------------------------------------------------------------------------------
#  PLT_DEBUG - Print plot data structure. 

procedure plt_debug (pl)

#  Arguments:
pointer	pl			# I: plot descriptor

begin
	call eprintf ("\n** Plot Data Structure Contents **\n\n")
	call eprintf (" Plot type: %d; diagnostic type: %d\n")
	    call pargi (PL_TYPE(pl))
	    call pargi (PL_DIAG_TYPE(pl))

	if (PL_GP(pl) == NULL)
	    call eprintf (" Plot graphics descriptor is NULL.\n")

	if (PL_AT(pl) == NULL)
	    call eprintf (" Plot atomic data descriptor is NULL.\n")

	call eprintf ("\nPlot Limits\n")
	call eprintf ("  Density: %8g:%8g\n  Temperature: %8.1f:%8.1f\n  Intensity: %8g:%8g\n")
	    call pargr (HI_NE(pl))
	    call pargr (LO_NE(pl))
	    call pargr (HI_TE(pl))
	    call pargr (LO_TE(pl))
	    call pargr (HI_INTENS(pl))
	    call pargr (LO_INTENS(pl))

	call eprintf ("\nPlot Contours\n")
	call eprintf ("  Primary curve: %8g; interval: %8g\n")
	    call pargr (PL_RATIO(pl))
	    call pargr (PL_DELTA(pl))
	call eprintf ("  List of secondary curves: %s\n")
	    call pargstr (PL_LIST(pl))
	call eprintf ("  No. contours: %d; No. points/curve: %d\n")
	    call pargi (PL_NCONTOUR(pl))
	    call pargi (PL_NPTS(pl))

	call eprintf ("  Color of primary: %d; secondary: %d\n")
	    call pargi (PL_REF_COLOR(pl))
	    call pargi (PL_DELTA_COLOR(pl))

	call eprintf ("\nText strings\n")
	call eprintf ("  User Title: %s\n  Plot Title: %s\n")
	    call pargstr (USER_TITLE(pl))
	    call pargstr (TITLE(pl))
	call eprintf ("  X-label: %s\n  Y-label: %s\n")
	    call pargstr (XLABEL(pl))
	    call pargstr (YLABEL(pl))
	call eprintf ("  Diagnostic Expression: %s\n  Output Table: %s\n")
	    call pargstr (DIAG_EXPR(pl))
	    call pargstr (PL_TABLE(pl))

	call flush (STDERR)
end

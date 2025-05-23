# Copyright restrictions apply - see st4gem$copyright.stsdas

include <mach.h>
include	"../atom.h" 
include	"../neberr.h" 
include	"ntplot.h" 

define	DEBUG	false

#---------------------------------------------------------------------9 Jun 97--
.help ntcl_init.x Jun97 nebular/plot
.ih
NAME
.nf
    get_cl_pltpars - Set up initial plot attributes from CL parameters
 get_basic_pltpars - Set up basic plot attributes from CL parameters
.fi
.endhelp
#-------------------------------------------------------------------------------
#  GET_CL_PLTPARS -- Set up initial plot attributes from CL parameters. 

procedure get_cl_pltpars (pl)

#  Calling arguments:
pointer	pl			# I: local plot descriptor

#  Declarations:
int	atom, ion		# atomic number, spectrum number
int	clgeti()		# fetch CL parameter TY_INT
real	clgetr()		# fetch CL parameter TY_REAL
int	clgwrd()		# fetch a keyword from an enumerated string
int	get_atom_typ()		# get atomic number from char string 
int	i			# generic
char	keywd[SZ_LINE]		# matched keyword
int	n_pts			# no. points calculated per curve
char	s1[SZ_LINE]		# generic
int	strlen()		# returns length of string
int	word_indx		# index of dictionary entry

errchk	at_alloc, diag_update, get_atom_typ, ratio_update

begin
	n_pts = clgeti ("resolution")
	if (n_pts < MIN_PTS) 
	    call error (BAD_PLOT_RESOL, "Insufficient resolution for plot")
	else
	    call plt_realloc (pl, n_pts)

	# Determine diagnostic type
	PL_TYPE(pl) = clgwrd ("plot_type", keywd, SZ_LINE, PLOT_TYPE)
	switch (PL_TYPE(pl)) {
	case INTENS_VS_NE:
	    PL_DIAG_TYPE(pl) = INTENSITY

	case INTENS_VS_TE:
	    PL_DIAG_TYPE(pl) = INTENSITY

	case TE_VS_NE:
	    PL_DIAG_TYPE(pl) = clgwrd ("diag_type", keywd, SZ_LINE, TEMDEN)

	default:
	    call error (BAD_NT, "Invalid plot type specified")
	}

	call get_basic_pltpars (pl, true)

	# Get Atomic number & ionization stage.
	# Note: ion number is zero-indexed, but task parameter is one-indexed.
	word_indx = clgwrd ("atom", keywd, SZ_LINE, ATOM)
	atom = get_atom_typ (keywd)
	ion = clgeti ("spectrum")
	ion = ion - 1

	# Allocate atomic data structure. 
	call at_alloc (atom, ion, PL_AT(pl))
	if (DEBUG)
	    call at_debug (PL_AT(pl))

	# Fetch transition description.
	call clgstr ("transition", USER_TRANSITION(pl), SZ_LINE)
	call diag_update (pl, USER_TRANSITION(pl), SZ_LINE)

	# Fetch primary curve value.
	call clgstr ("primary", s1, SZ_LINE)
	call ratio_update (pl, s1)

	# Populate contour list; construct one if none is given.
	call clgstr ("clist", PL_LIST(pl), SZ_LINE)
	PL_NCONTOUR(pl) = clgeti ("ncontours")
	PL_DELTA(pl)    = clgetr ("interval")

	if (strlen (PL_LIST(pl)) == 0 && PL_NCONTOUR(pl) > 0) {
	    if (PL_DELTA(pl) > EPSILONR) {
	    	for (i = 1; i <= PL_NCONTOUR(pl); i = i+1) {
		    call sprintf (s1, SZ_LINE, "%g %g ")
	    	    	call pargr (PL_RATIO(pl) * 10**(  i * PL_DELTA(pl)))
	    	    	call pargr (PL_RATIO(pl) * 10**(-(i * PL_DELTA(pl))))
		    call strcat (s1, PL_LIST(pl), SZ_LINE)
		}
	    }
	}

	# Fetch user title & curve colors.
	call clgstr ("title", USER_TITLE(pl), SZ_LINE)
	PL_REF_COLOR(pl) = clgeti ("ref_color")
	PL_DELTA_COLOR(pl) = clgeti ("delta_color")
end

 
#-------------------------------------------------------------------------------
#  GET_BASIC_PLTPARS -- Set up basic plot attributes from CL parameters. 

procedure get_basic_pltpars (pl, get_intens)

#  Calling arguments:
pointer	pl			# I: local plot descriptor
bool	get_intens		# I: get intensity limits?

#  Declarations:
bool	clgetb()		# fetch CL parameter TY_BOOL
real	clgetr()		# fetch CL parameter TY_REAL

errchk	

begin
	# Set initial plot limits
	LO_NE(pl) = clgetr ("min_dens")
	HI_NE(pl) = clgetr ("max_dens")
	if (LO_NE(pl) >= HI_NE(pl))
	    call error (BAD_NT, "Minimum N_e must exceed maximum!")
	PLOT_LOG_NE(pl) = clgetb ("log_ne")

	LO_TE(pl) = clgetr ("min_temp")
	HI_TE(pl) = clgetr ("max_temp")
	if (LO_TE(pl) >= HI_TE(pl))
	    call error (BAD_NT, "Minimum T_e must exceed maximum!")
	PLOT_LOG_TE(pl) = clgetb ("log_te")

	LO_INTENS(pl) = INDEFR
	HI_INTENS(pl) = INDEFR
	if (get_intens) {
	    LO_INTENS(pl) = clgetr ("min_intens")
	    HI_INTENS(pl) = clgetr ("max_intens")
	    if (LO_INTENS(pl) >= HI_INTENS(pl))
	    	call error (BAD_NT, "Minimum intensity must exceed maximum!")
	    PLOT_LOG_INTENS(pl) = clgetb ("log_intens")

	} 
end

 

# Copyright restrictions apply - see st4gem$copyright.stsdas

include	"../at.h"
include	"../atom.h"
include	"../neberr.h"
include	"ntplot.h" 

define	CMDS	"|atom|clist|diagtype|dmin|dmax|dscale|imin|imax|iscale|\
	|interval|ion|ncontours|plottype|primary|refcolor|resolution|title|tmin|\
		|tmax|transition|tscale|write"
define	SZ_CMD		31
define	Command		Memc[cmd]
define	LOG_LIN		"|log|linear|"

define	SET_ATOM	1
define	CLIST		2
define	DIAGTYPE	3
define	DMIN		4
define	DMAX		5
define	DSCALE		6
define	IMIN		7
define	IMAX		8
define	ISCALE		9
# skip			10
define	INTERVAL	11
define	SET_ION		12
define	NCONTOURS	13
define	PLOTTYPE	14
define	PRIMARY		15
define	REFCOLOR	16
define	RESOLUTION	17
define	TITLE		18
define	TMIN		19
# skip			20
define	TMAX		21
define	TRANSITION	22
define	TSCALE		23
define	WRITE		24

#---------------------------------------------------------------------6 Jul 97--
.help nt_colon.x Apr97 nebular/plot
.ih
NAME
.nf
  nt_colon - Interpret colon commands
.fi
.endhelp
#-------------------------------------------------------------------------------
#  NT_COLON -- Interpret colon commands.

procedure nt_colon (pl, o)

#  Calling arguments:
pointer	pl			# I: plot state descriptor
pointer	o			# I: parameter object

#  Declarations:
pointer	at			# atomic data object
bool	at_flag			# was atom or ion updated?
pointer	cmd			# colon command
double	dval			# generic
int	get_atom_typ()		# get atomic number from atom name
int	i, ival, ival2		# generic
int	ncmd			# index of command
bool	pr_npname()		# get name of next parameter
double	pr_getd()		# get parameter value, TY_DOUBLE
int	pr_geti()		# get parameter value, TY_INT
char	s1[SZ_LINE]		# generic
int	strdic()		# search a string dictionary
bool	streq()			# are two strings equal?
int	strlen()		# returns length of string

errchk	at_alloc, diag_update, pr_getd, pr_geti, pr_npname

begin
	at_flag = false
	at = PL_AT(pl)

	# Process each specified parameter.
	while (pr_npname (o, cmd)) {
	    ncmd = strdic (Command, s1, SZ_CMD, CMDS)
	    switch (ncmd) {

	    # Select a different atom:
	    case SET_ATOM:
	    	call pr_gets (o, Command, s1, SZ_CMD)
	    	ival = strdic (s1, s1, SZ_CMD, ATOM)
		ival = get_atom_typ (s1)
	    	if (ival > 0) {
		    at_flag = true
		    AT_ATOM(at) = ival

		} else 
	    	    call printf ("Invalid atom specified\7\n")

	    # Select a different ion:
	    case SET_ION:
	    	ival = pr_geti (o, Command)
	    	if (ival >= 1) {
		    at_flag = true
		    AT_ION(at) = ival - 1

		} else 
	    	    call printf ("Invalid ion specified\7\n")

	    # Set secondary contour list:
	    case CLIST:
	    	call pr_gets (o, Command, PL_LIST(pl), SZ_CMD)
	    	if (strlen (PL_LIST(pl)) == 0 && PL_NCONTOUR(pl) > 0) {
		    for (i = 1; i <= PL_NCONTOUR(pl); i = i+1) {
			call sprintf (s1, SZ_LINE, "%g %g ")
		    	    call pargr (PL_RATIO(pl) * 10**(  i * PL_DELTA(pl)))
		    	    call pargr (PL_RATIO(pl) * 10**(-(i * PL_DELTA(pl))))
			call strcat (s1, PL_LIST(pl), SZ_LINE)
		    }
		} 

	    # Change the diagnostic type:
	    case DIAGTYPE:
	    	call pr_gets (o, Command, s1, SZ_CMD)
		call strlwr (s1)
	    	ival = strdic (s1, s1, SZ_CMD, TEMDEN)
	    	if (ival > 0 && PL_TYPE(pl) == TE_VS_NE) 
		    PL_DIAG_TYPE(pl) = ival
		else 
	    	    call printf ("Invalid diagnostic type specified\7\n")

	    # Change the density limits:
	    case DMIN:
	    	dval = pr_getd (o, Command)
	    	if (dval >= 1.D+0 && dval <= 1.d+8)
		    LO_NE(pl) = dval

	    case DMAX:
	    	dval = pr_getd (o, Command)
	    	if (dval >= 1.D+0 && dval <= 1.d+8)
		    HI_NE(pl) = dval

	    # Select log/linear plotting for density:
	    case DSCALE:
	    	call pr_gets (o, Command, s1, SZ_CMD)
	    	ival = strdic (s1, s1, SZ_CMD, LOG_LIN)
	    	if (streq (s1, "log") )
		    PLOT_LOG_NE(pl) = true

	    	else if (streq (s1, "linear") )
		    PLOT_LOG_NE(pl) = false

		else 
	    	    call printf ("Invalid log/linear scale specified\7\n")

	    # Change the intensity limits:
	    case IMIN:
	    	dval = pr_getd (o, Command)
	    	if (dval > 1.D-7)
		    LO_INTENS(pl) = dval

	    case IMAX:
	    	dval = pr_getd (o, Command)
	    	if (dval >= 1.D+0)
		    HI_INTENS(pl) = dval

	    # Change the interval between contours:
	    case INTERVAL:
	    	dval = pr_getd (o, Command)
	    	if (dval > 1.D-7)
		    PL_DELTA(pl) = dval

	    # Select log/linear plotting for intensity:
	    case ISCALE:
	    	call pr_gets (o, Command, s1, SZ_CMD)
	    	ival = strdic (s1, s1, SZ_CMD, LOG_LIN)
	    	if (streq (s1, "log") )
		    PLOT_LOG_INTENS(pl) = true

	    	else if (streq (s1, "linear") )
		    PLOT_LOG_INTENS(pl) = false

		else 
	    	    call printf ("Invalid log/linear scale specified\7\n")

	    # Change the number of contours:
	    case NCONTOURS:
	    	ival = pr_geti (o, Command)
	    	if (ival >= 0)
		    PL_NCONTOUR(pl) = ival

		else {
	    	    call printf ("Invalid no. contours specified: %d\7\n")
		    call pargi (ival)
		}

	    # Change the plot type:
	    case PLOTTYPE:
	    	call pr_gets (o, Command, s1, SZ_CMD)
		call strupr (s1)
	    	ival = strdic (s1, s1, SZ_CMD, PLOT_TYPE)
	    	if (ival > 0 && ival < 4) {
		    PL_TYPE(pl) = ival
		    if (ival == INTENS_VS_NE || ival == INTENS_VS_TE)
		    	PL_DIAG_TYPE(pl) = INTENSITY

		    else {
			if (PL_DIAG_TYPE(pl) != DENSITY && 
			    PL_DIAG_TYPE(pl) != TEMPERATURE)
		    	    PL_DIAG_TYPE(pl) = DENSITY
		    }

		} else {
	    	    call printf ("Invalid plot type specified: %s = %d\7\n")
		    call pargstr (s1)
		    call pargi (ival)
		}

		call strcpy ("default", s1, SZ_LINE)
		call diag_update (pl, s1, SZ_LINE)

	    # Change expression for primary curve value:
	    case PRIMARY:
	    	call pr_gets (o, Command, s1, SZ_LINE)
		call ratio_update (pl, s1)

	    # Change the color of the reference curve:
	    case REFCOLOR:
	    	ival = pr_geti (o, Command)
	    	if (ival >= 0 && ival < 16)
		    PL_REF_COLOR(pl) = ival

		else {
	    	    call printf ("Invalid reference color specified: %d\7\n")
		    call pargi (ival)
		}

	    # Change the plot resolution:
	    case RESOLUTION:
	    	ival = pr_geti (o, Command)
	    	if (ival >= MIN_PTS) 
		    call plt_realloc (pl, ival)

		else {
	    	    call printf ("Invalid resolution specified: %d\7\n")
		    call pargi (ival)
		}

	    # Set plot title:
	    case TITLE:
	    	call pr_gets (o, Command, s1, SZ_LINE)
	    	if (streq (s1, "default") )
	    	    call strcpy (EOS, USER_TITLE(pl), SZ_LINE)
		else 
	    	    call strcpy (s1, USER_TITLE(pl), SZ_LINE)

		call title_update (pl)

	    # Change the temperature limits:
	    case TMIN:
	    	dval = pr_getd (o, Command)
	    	if (dval >= 1.D+0 && dval <= 1.d+8)
		    LO_TE(pl) = dval

	    case TMAX:
	    	dval = pr_getd (o, Command)
	    	if (dval >= 1.D+0 && dval <= 1.d+8)
		    HI_TE(pl) = dval

	    # Set custom transition:
	    case TRANSITION:
	    	call pr_gets (o, Command, s1, SZ_LINE)
		call diag_update (pl, s1, SZ_LINE)

	    # Select log/linear plotting for temperature:
	    case TSCALE:
	    	call pr_gets (o, Command, s1, SZ_CMD)
	    	ival = strdic (s1, s1, SZ_CMD, LOG_LIN)
	    	if (streq (s1, "log") )
		    PLOT_LOG_TE(pl) = true

	    	else if (streq (s1, "linear") )
		    PLOT_LOG_TE(pl) = false

		else 
	    	    call printf ("Invalid log/linear scale specified\7\n")

	    # Write the curves to an output table:
	    case WRITE:
	    	call printf ("Curve output not yet supported.\7\n")
#	    	call pr_gets (o, Command, s1, SZ_CMD)
#		tp = tbtopn (s1, READ_WRITE, NULL)
#		call curve_to_table (tp, c_name, diagtype, array, npts)
#		call tbtclo (tp)

	    # Unrecognized cursor command -- ring bell.
	    default:
	    	call printf ("Invalid colon command: %s\7\n")
		call pargstr (Command)
	    }
	}

	# Reset atom/ion if needed
	if (at_flag) {
	    ival  = AT_ATOM(at)
	    ival2 = AT_ION(at)
	    call at_alloc (ival, ival2, at)
	    call strcpy ("default", s1, SZ_LINE)
	    call diag_update (pl, s1, SZ_LINE)
	}
end



# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<tbset.h>
include	"../atom.h" 
include	"../neberr.h" 

define	ATOM_COL	1
define	ION_COL		2
define	ZONE_COL	3
define	NAME_COL	4
define	DIAG_COL	5

define	N_TCOLS		5

#--------------------------------------------------------------------14 Jul 97--
.help curve_ref.x Mar97 nebular/plot
.ih
NAME
.nf
   open_curve - Open the diagnostic plot data table
get_diag_line - Fetch diagnostic curve info from table row
.fi
.endhelp
#-------------------------------------------------------------------------------
#  OPEN_CURVE -- Open the diagnostic curve data table.

procedure open_curve (tp, colptr, n_rows)

#  Calling arguments:
pointer	tp			# O: table descriptor
pointer	colptr[N_TCOLS]		# O: column pointers
int	n_rows			# O: number of table rows

#  Declarations:
char	colname[SZ_COLNAME]	# table column names
char	errmsg[SZ_LINE]		# error message
int	i			# generic
char	path[SZ_FNAME]		# path to atomic data files
int	tbpsta()		# fetch table attributes

string	table_name  "ntplot.dat"

begin
	# Open table containing the emission line descriptions
	call open_catalog (table_name, tp, path)

	# Find columns by name, or report an error. 
	do i = 1, N_TCOLS {
	    call sprintf (colname, SZ_COLNAME, "c%1d")
		call pargi (i)
	    call tbcfnd (tp, colname, colptr[i], 1)

	    if (colptr[i] == NULL) {
		call sprintf (errmsg, SZ_LINE, 
			"Unable to find column %s in reference table")
		call pargstr (colname)
		call error (BAD_REF_DATA, errmsg)
	    }
	}

	n_rows = tbpsta (tp, TBL_NROWS)
end


#-------------------------------------------------------------------------------
#  GET_DIAG_LINE - Fetch diagnostic curve info from table row. 

procedure get_diag_line (tp, cp, row, atom, ion, zone, diag_type, colname)

#  Arguments:
pointer	tp			# I: table descriptor
pointer	cp[N_TCOLS]		# I: column descriptor array
int	row			# I: table row number
int	atom			# O: atomic number 
int	ion			# O: spectrum (ion) number
int	zone			# O: ionization zone
int	diag_type		# O: diagnostic type
char	colname[ARB]		# O: flux column name 

# Declarations:
int	get_atom_typ()		# match atom name to symbolic 
char	s1[SZ_COLNAME]		# generic 
int	strdic()		# search a string dictionary

begin
	# Fetch each element
	call tbegtt (tp, cp[ATOM_COL], row, s1, SZ_COLNAME)
	atom = get_atom_typ (s1)

	call tbegtt (tp, cp[DIAG_COL],  row, s1, SZ_COLNAME)
	call strlwr (s1)
	diag_type = strdic (s1, s1, SZ_COLNAME, TEMDEN)

	call tbegti (tp, cp[ION_COL],   row, ion)
	call tbegti (tp, cp[ZONE_COL],  row, zone)
	call tbegtt (tp, cp[NAME_COL],  row, colname, SZ_COLNAME)

end



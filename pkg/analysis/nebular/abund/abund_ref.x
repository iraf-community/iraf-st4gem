# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<tbset.h>
include	"../neberr.h" 

define	ATOM_COL	1
define	ION_COL		2
define	ZONE_COL	3
define	WAVE_COL	4
define	WT_COL		5
define	NLINE_COL	6
define	WIDTH_COL	7
define	NAME_COL	8

#--------------------------------------------------------------------14 Jul 97--
.help abund_ref.x Feb97 nebular/fivel
.ih
NAME
.nf
    open_abund - Open the abundance data table
get_abund_line - Fetch a row of abundance line data
.fi
.endhelp
#-------------------------------------------------------------------------------
#  OPEN_ABUND -- Open the abundance data table.

procedure open_abund (tp, colptr, n_rows)

#  Calling arguments:
pointer	tp			# O: table descriptor
pointer	colptr[8]		# O: column pointers
int	n_rows			# O: number of table rows

#  Declarations:
char	colname[SZ_COLNAME]	# table column names
char	errmsg[SZ_LINE]		# error message
int	i			# generic
char	path[SZ_FNAME]		# path to atomic data files
int	tbpsta()		# fetch table attributes

string	table_name  "abund.dat"

begin
	# Open table containing the emission line descriptions
	call open_catalog (table_name, tp, path)

	# Find columns by name, or report an error. 
	do i = 1, 8 {
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
#  GET_ABUND_LINE - Fetch a row of abundance line data. 

procedure get_abund_line (tp, cp, row, atom, ion, zone, wave, wt, 
				nlines, width, colname)

#  Arguments:
pointer	tp		# I: table descriptor
pointer	cp[8]		# I: column descriptor array
int	row		# I: table row number
int	atom		# O: atomic number 
int	ion		# O: spectrum (ion) number
int	zone		# O: ionization zone
real	wave		# O: wavelength
real	wt		# O: weight
int	nlines		# O: no. lines in multiplet
real	width		# O: wavelength interval of multiplet
char	colname[ARB]	# O: flux column name 

begin

	# Fetch each element
	call tbegti (tp, cp[ATOM_COL],  row, atom)
	call tbegti (tp, cp[ION_COL],   row, ion)
	call tbegti (tp, cp[ZONE_COL],  row, zone)
	call tbegtr (tp, cp[WAVE_COL],  row, wave)
	call tbegtr (tp, cp[WT_COL],    row, wt)
	call tbegti (tp, cp[NLINE_COL], row, nlines)
	call tbegtr (tp, cp[WIDTH_COL], row, width)
	call tbegtt (tp, cp[NAME_COL],  row, colname, SZ_COLNAME)

end



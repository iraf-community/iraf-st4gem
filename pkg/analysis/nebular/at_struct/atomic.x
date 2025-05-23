# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<error.h>
include	<tbset.h>
include	"../at.h"
include	"../fivel.h"
include	"../neberr.h"

define	ELVL_COL	1
define	WT_COL		2
define	RAD_COL		3

#--------------------------------------------------------------------14 Jul 97--
.help atomic.x Jan96 nebular/fivel
.ih
NAME
.nf
     read_atomic - Read the T_e - independent atomic data from an open table
     at_col_find - Fetch column pointers by names for atomic data
  open_atom_file - Open the appropriate atomic datafile for the specified ion
calc_energy_diff - Compute atomic energy level differences 
.fi
.endhelp
#-------------------------------------------------------------------------------
#  READ_ATOMIC -- Read the T_e - independent atomic data from an open table. 

procedure read_atomic (tp, at)

#  Calling arguments:
pointer	tp		# I: table descriptor
pointer	at		# I: atomic data object

#  Local variables:
int	col		# index of table column, row
pointer	colptr		# table column pointers
int	last_col	# last table column to read
int	nlvl		# no. energy levels for this ion
int	ncols		# no. table columns
pointer	nullflag	# array of null flags
pointer	sp		# top of stack

string	col_template	"A(*->%d)"	# table column name template

#  Functions used:
int	tbpsta()	# fetch table attributes

#  Memory management:
define	Cptr		Memi[colptr+$1-1]
define	Nulflg		Memb[nullflag]

define	RAD_TR_PROB	Memd[AT_RAD_TR(at)+($1-1)+(($2-1)*AT_NLVL(at))]

errchk	at_col_find

begin
	n_lvl = AT_NLVL(at)
	ncols = tbpsta (tp, TBL_NCOLS)

	call smark (sp)
	call salloc (colptr, ncols, TY_POINTER)
	call salloc (nullflag, n_lvl, TY_BOOL)

	# Initialize column names. 
	call at_col_find (tp, col_template, Cptr(1), n_lvl)

	# Get energy levels & statistical weights from the table.
	call tbcgtd (tp, Cptr(ELVL_COL), E_TRANS(at), Nulflg, 1, n_lvl)
	call tbcgti (tp, Cptr(WT_COL),   WEIGHT(at),  Nulflg, 1, n_lvl)

	# Populate the radiative transition probabilities array.
	# Note that only the columns through A(*->4) need be read.
	last_col = n_lvl + RAD_COL - 2
	do col = RAD_COL, last_col
	    call tbcgtd (tp, Cptr(col), RAD_TR_PROB(1,col-2), Nulflg, 1, n_lvl)

	call sfree (sp)
end


#-------------------------------------------------------------------------------
#  AT_COL_FIND -- Fetch column pointers by names for atomic data.  

procedure at_col_find (tp, template, colptr, nlvls)

#  Calling arguments:
pointer	tp			# I: table descriptor
char	template[ARB]		# I: table column name template
pointer	colptr[ARB]		# O: table column pointers
int	nlvls			# I: no. atomic energy levels expected

#  Local variables:
char	errmsg[SZ_LINE]		# error message
int	i, j			# generic
char	name[SZ_COLNAME]	# 

begin
	# Find transition energy column.
	call tbcfnd (tp, "Energy", colptr[ELVL_COL], 1)
	if (colptr[ELVL_COL] == NULL) {
	    call sprintf (errmsg, SZ_LINE, 
		"Unable to find transition energy column in reference table")
	    call error (BAD_REF_DATA, errmsg)
	}

	# Find statistical weight column.
	call tbcfnd (tp, "Stat_Weight", colptr[WT_COL], 1)
	if (colptr[WT_COL] == NULL) {
	    call sprintf (errmsg, SZ_LINE, 
		"Unable to find statistical weight column in reference table")
	    call error (BAD_REF_DATA, errmsg)
	}

	# Find radiative transition prob. columns by name. 
	do i = 1, nlvls-1 {
	    j = i + RAD_COL - 1
	    call sprintf (name, SZ_COLNAME, template)
		call pargi (i)
	    call tbcfnd (tp, name, colptr[j], 1)
	    if (colptr[j] == NULL) {
		call sprintf (errmsg, SZ_LINE, 
			"Unable to find column %s in reference table")
		call pargstr (name)
		call error (BAD_REF_DATA, errmsg)
	    }
	}

end


#-------------------------------------------------------------------------------
#  OPEN_ATOM_FILE -- Open the appropriate atomic data table for the specified 
#			ion.

pointer procedure open_atom_file (atom, ion, at_data_type)

#  Calling arguments:
int	atom			# I: atomic number
int	ion			# I: spectrum number
int	at_data_type		# I: atomic data file type
pointer	tp			# O: table descriptor

#  Declarations:
int	at_match_row()		# find matching row in input tabl
char	colname[SZ_COLNAME]	# table column names
pointer	colptr[4]		# column pointers
char	errmsg[SZ_LINE]		# error message
int	get_atom_str()		# get atom/ion symbol
int	i			# generic
char	path[SZ_FNAME]		# path to atomic data files
int	row			# table row
char	text[SZ_FNAME]		# generic
char	tbl_name[SZ_FNAME]	# atomic data file to open
pointer	tbtopn()		# open named table

errchk	open_catalog

begin
	# Open table containing the names of the files
	call open_catalog ("at_refer.dat", tp, path) 

	# Find columns by name, or report an error. 
	do i = 1, 4 {
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

	# Fetch the appropriate file name.
	row = at_match_row (tp, colptr[1], colptr[2], atom, ion)
	if (row < 1) {
	    i = get_atom_str (atom, ion, text, SZ_FNAME)
	    call sprintf (errmsg, SZ_LINE, 
		"Unable to locate file in atomic reference table for atom %s")
		call pargstr (text)

	    call error (REF_TABLE_ACCESS, errmsg)
	}

	call tbegtt (tp, colptr[at_data_type+2], row, text, SZ_FNAME)
	call tbtclo (tp)

	# Open the table & return the pointer
	call strcpy (path, tbl_name, SZ_FNAME)
	call strcat (text, tbl_name, SZ_FNAME)
	iferr (tp = tbtopn (tbl_name, READ_ONLY, NULL) ) {
	    call sprintf (errmsg, SZ_LINE, 
		"Unable to open atomic data table: %s")
		call pargstr (tbl_name)
	    call error (REF_TABLE_ACCESS, errmsg)
	} else
	    return (tp)

end


#-------------------------------------------------------------------------------
#  AT_MATCH_ROW - Find the row matching the atomic & spectrum numbers in the 
#		input table. 

int procedure at_match_row (tp, cp_atom, cp_ion, atom, ion)

#  Arguments:
pointer	tp		# I: table descriptor
pointer	cp_atom		# I: column descriptor for atom
pointer	cp_ion		# I: column descriptor for ion
int	atom		# I: atomic number 
int	ion		# I: spectrum (ion) number
int	row		# O: selected row

#  Local variables:
int	val1, val2	# generic
int	n_rows		# number of rows in table
int	tbpsta()	# table status

define	NO_MATCH	0

begin
	n_rows = tbpsta (tp, TBL_NROWS)

	# Match to first successful atomic/spectrum number pair
	do row = 1, n_rows {
	    call tbegti (tp, cp_atom, row, val1)
	    call tbegti (tp, cp_ion, row, val2)
	    if (atom == val1 && ion == val2)
	    	return (row)
	}

	return (NO_MATCH)
end


#-------------------------------------------------------------------------------
#  CALC_ENERGY_DIFF -	Compute atomic energy level differences given energy of 
#			each level above ground state.  

procedure calc_energy_diff (energy_lvl, energy_diff, n_lvl)

#  Arguments
double	energy_lvl[n_lvl]		# I: energy above ground state
double	energy_diff[n_lvl,n_lvl]	# O: relative energy separation
int	n_lvl				# I: no. atomic energy levels

#  Local variables
int	i, j				# generic

define	HC	1.986485d-8	# product of Planck's constant & speed of light

begin
	do i = 1, n_lvl {
	    do j = 1, i
		energy_diff[j,i] = 0.

	    do j = i+1, n_lvl
		energy_diff[j,i] = HC * (energy_lvl[j] - energy_lvl[i])
	}

end



# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<tbset.h>
include	<mach.h>
include	"../fivel.h"
include	"../neberr.h"

define	MIN_REAL	1./MAX_REAL

#--------------------------------------------------------------------11 Jul 97--
.help tab_io.x Jan96 nebular/fivel
.ih
NAME
  open_catalog - Open a catalog in the atomic data reference directory. 
curve_to_table - Write diagnostic curve to ouput table
      get_flux - Fetch flux values from the input table. 
get_extinction - Fetch extinction law & constant input table. 
.endhelp
#-------------------------------------------------------------------------------
#  OPEN_CATALOG - Open a catalog in the atomic data reference directory. 

procedure open_catalog (catalog, tp, path)

#  Calling arguments
char	catalog[ARB]		# I: catalog table to open
pointer	tp			# O: returned value
char	path[SZ_FNAME]		# O: path to atomic data files

#  Local variables
int	envgets()		# fetch IRAF environment variable TY_STRING
int	nchar			# dummy
char	tbl_name[SZ_FNAME]	# atomic data file to open
pointer	tbtopn()		# open named table

begin
	# Construct the pathname to the atomic data files.
	call clgstr ("at_data", path, SZ_FNAME)
	nchar = envgets (path, path, SZ_FNAME)

	# Open table containing the names of the files
	call strcpy (path, tbl_name, SZ_FNAME)
	call strcat (catalog, tbl_name, SZ_FNAME)
	iferr (tp = tbtopn (tbl_name, READ_ONLY, NULL) ) {
	    call eprintf ("Error reading atomic data file: %s\n")
	    	call pargstr (tbl_name)
	    call flush (STDERR)
	    call error (REF_TABLE_ACCESS, "Error reading atomic reference data")
	}
end


#-------------------------------------------------------------------------------
#  CURVE_TO_TABLE - Write diagnostic curve to ouput table. 

procedure curve_to_table (tp, c_name, diag_type, array, npts)

#  Arguments:
pointer	tp			# I: output table descriptor
char	c_name[SZ_COLNAME]	# I: output column name
int	diag_type		# I: diagnostic type (DENSITY|TEMPERATURE)
real	array[ARB]		# I: array to write
int	npts			# I: size of array

#  Declarations:
char	c_fmt[SZ_COLFMT]	# output column format
pointer	c_ptr			# output column pointer
char	c_unit[SZ_COLUNITS]	# output column units

include	"../atom.h"

begin
	if (diag_type == DENSITY) {
	    call strcpy ("1/cm^3", cunit, SZ_COLUNITS)
	    call strcpy ("%10.1f", cfmt, SZ_COLFMT)

	} else if (diag_type == TEMPERATURE) {
	    call strcpy ("Kelvin", cunit, SZ_COLUNITS)
	    call strcpy ("%8.0f", cfmt, SZ_COLFMT)

	} else if (diag_type == INTENSITY) {
	    call strcpy ("none", cunit, SZ_COLUNITS)
	    call strcpy ("%8.3f", cfmt, SZ_COLFMT)
	}

	if (tp != NULL) {
	    call tbcdef (tp, c_ptr, c_name, c_unit, c_fmt, TY_REAL, 1, 1)
	    call tbcptr (tp, c_ptr, array, 1, n_pts)
	}

end


#-------------------------------------------------------------------------------
#  GET_FLUX -	Fetch values from the input table.  If the line flux cannot 
#		be found or is <= zero, return INDEF. 

real procedure get_flux (itp, row, c_name)

#  Calling arguments
pointer	itp			# I: input table descriptor
int	row			# I: input table row
char	c_name[ARB]		# I: user column name
real	val			# O: returned value

#  Local variables
char	colname[SZ_COLNAME+1]	# name of table column
pointer	colptr			# column descriptor for input table

begin
	call clgstr (c_name, colname, SZ_COLNAME)
	call tbcfnd (itp, colname, colptr, 1)

	if (colptr != NULL) {
	    call tbegtr (itp, colptr, row, val)

	} else 
	    val = INDEFR

	if (val < MIN_REAL)
	    val = INDEFR

	return (val)
end


#-------------------------------------------------------------------------------
#  GET_EXTINCTION - Fetch extinction law & constant input table. 

procedure get_extinction (itp, row, ext_fn, c_ext)

#  Arguments:
pointer	itp			# I: input table descriptor
int	row			# I: input table row
int	ext_fn			# O: index of I.S. extinction function
real	c_ext			# O: extinction constant

#  Local variables:
char	colname[SZ_COLNAME+1]	# name of column
pointer	colptr			# column descriptor for input table
bool	ext_flag		# are input fluxes extinction corrected?
char	red_law[3]		# choice of reddening law

#  Functions used:
int	strdic()		# choose a string from a dictionary

begin
	# Fetch the extinction correction flag. 

	call clgstr ("extcorr_col", colname, SZ_COLNAME)
	call tbcfnd (itp, colname, colptr, 1)
	if (colptr == NULL) 
	    ext_flag = true
	else 
	    call tbegtb (itp, colptr, row, ext_flag)

	# Fetch the extinction constant; default value is 0. 

	call clgstr ("c_ext_col", colname, SZ_COLNAME)
	call tbcfnd (itp, colname, colptr, 1)
	if (colptr == NULL) 
	    c_ext = 0.0
	else 
	    call tbegtr (itp, colptr, row, c_ext)

	if (c_ext < MIN_REAL || ext_flag || IS_INDEF (c_ext))
	    c_ext = 0.0

	# Default to the GALACTIC extinction law unless another legitimate 
	# choice is given. 
	call clgstr ("redlaw_col", colname, SZ_COLNAME)
	call tbcfnd (itp, colname, colptr, 1)

	if (colptr == NULL) 
	    ext_fn = GAL
	else {
	    call tbegtt (itp, colptr, row, red_law, 3)
	    call strlwr (red_law)
	    ext_fn = strdic (red_law, red_law, 3, EXTN_MODEL)
	    if (ext_fn <= 0)
	    	ext_fn = GAL
	}

end




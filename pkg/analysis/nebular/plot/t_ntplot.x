# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<error.h>
include	<tbset.h>
include	"../fivel.h"
include	"../neberr.h"

define	DEBUG	false

#---------------------------------------------------------------------7 Sep 97--
.help t_ntplot.x Feb97 nebular/fivel
.ih
NAME
  t_ntplot - Task to plot density vs. temperature diagnostic curves. 
.endhelp
#-------------------------------------------------------------------------------
#  T_NTPLOT -	Task to plot N_e vs. T_e diagnostic diagram. 

procedure t_ntplot ()

#  Declarations:
pointer	obj_col, reg_col	# column descriptors from input table
char	errmsg[SZ_LINE] 	# text of error message
char	fluxtab[SZ_FNAME] 	# name for table of line fluxes
pointer	ilr			# pointer to line ratio structure
pointer	itp, otp		# input, output table descriptors
int	mat_object()		# find matching object in flux table
char	object[SZ_FNAME]	# list of object names
char	outtab[SZ_FNAME] 	# table of output curves
char	region[SZ_FNAME]	# Region name from input table
int	row			# input table row 
bool	strne()			# are two strings not equal?
pointer	tbtopn()		# initialize table structure

begin
#	if (DEBUG)
#	    call memlog ("---------- Start ntplot ----------")

	# Get the input table name & plot limits
	call clgstr ("fluxtab",   fluxtab, SZ_FNAME)
	call imgcluster (fluxtab, fluxtab, SZ_FNAME)

	# Open the input table of line fluxes. 
	iferr (itp = tbtopn (fluxtab, READ_WRITE, NULL)) {
	    call sprintf (errmsg, SZ_LINE, "Input table %s cannot be accessed")
		call pargstr (fluxtab)
	    call error (TABLE_ACCESS, errmsg)

	# If successful, locate the Object_ID and Region columns.
	} else {
	    call tbcfnd (itp, "Object_ID", obj_col, 1)
	    call tbcfnd (itp, "Region",    reg_col, 1)
	}

	# Open output table for results if a name is specified. 
	call clgstr ("outtab", outtab, SZ_FNAME)
	if (strne (outtab, EOS)) {
	    otp = tbtopn (outtab, NEW_FILE, NULL)
	    call tbtcre (otp)
	} else
	    otp = NULL

	# Search the table to find matching row number. 
	call clgstr ("object", object, SZ_FNAME)
	call clgstr ("region", region, SZ_FNAME)
	row = mat_object (itp, obj_col, reg_col, object, region, false)
	if (row <= 0)
	    call error (BAD_OBJ_ID, "No matching Object_ID in input table")

	# Create data structure for input line ratios. 
	call malloc (ilr, LEN_ILR, TY_STRUCT)

	# Calculate line ratios and Plot diagnostic curves for object selected. 
	iferr {
	    call l_ratios (itp, ilr, row)
	    call diagplot (ilr, otp, object)

	} then {
	    call mfree (ilr, TY_STRUCT)
	    call tbtclo (itp)
	    call erract (EA_WARN)
	    call error (1, "Error generating diagnostic plot")
	}

	# Close input table and free memory. 
	call tbtclo (itp)
	if (otp != NULL)
	    call tbtclo (otp)
	call mfree (ilr, TY_STRUCT)

#	if (DEBUG)
#	    call memlog ("---------- End ntplot ----------")
end



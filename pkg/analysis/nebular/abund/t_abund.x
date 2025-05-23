# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<tbset.h>
include	<error.h>
include	"../fivel.h"
include	"../neberr.h"
include	"../neb_io.h"
include	"../zones.h"

#---------------------------------------------------------------------9 Jul 97--
.help t_abund.x Jul97 nebular/abund
.ih
NAME
.nf
    t_abund - Driver for ionic abundances program 
    find_tn - Fetch column ptrs for T_e/N_e diagnostics for @zone 
get_zone_nt - Fetch T_e & N_e diagnostics for each zone from input table
.fi
.endhelp
#-------------------------------------------------------------------------------
#  T_ABUND -	SPP driver for ionic abundances program.  Performs table 
#		I/O and fetches task parameters.  

procedure t_abund ()

#  Declarations
char	colname[SZ_COLNAME]	# Object/Region column names
pointer	did_ptr			# ID column pointer from diag table
pointer	dreg_ptr		# Region column pointer from diag table
real	dens[N_ZONES]		# densities for each zone
char	diagtab[SZ_FNAME] 	# name for table of T_e & N_e
int	drow, frow		# diagtab, fluxtab table row 
pointer	dtp, ftp		# diag, flux table descriptors
char	errmsg[SZ_LINE]		# error message
pointer	fid_ptr			# ID column pointer from flux table
pointer	freg_ptr		# Region column pointer from flux table
char	fluxtab[SZ_FNAME] 	# name for table of line fluxes
pointer	lfl			# pointer to flux object list
pointer	lfl_alloc()		# allocate flux object list
int	match_object()		# get flux table row that matches diag row
int	n_rows			# no. diag table rows
char	object[SZ_FNAME]	# Object_ID, diag file
bool	one_intable		# are the input & output files the same?
char	region[SZ_FNAME]	# Region code, diag file
bool	streq()			# are two strings equal?
int	tbpsta()		# get table header parameter
pointer	tbtopn()		# initialize table structure
real	temp[N_ZONES]		# temperatures for each zone

errchk	abund, getrow

begin
	# Get the input/output table names
	call clgstr ("intable", fluxtab, SZ_FNAME)
	call clgstr ("diagtab", diagtab, SZ_FNAME)

	# Open the input table of line fluxes; get table dimensions
	iferr (ftp = tbtopn (fluxtab, READ_WRITE, NULL)) 
	    call error (TABLE_ACCESS, "Error loading flux table")

	# Check to see if flux table = diag table
	one_intable = streq (fluxtab, diagtab)

	# If not, open the output table
	if (one_intable) 
	    dtp = ftp
	else {
	    iferr (dtp = tbtopn (diagtab, READ_WRITE, NULL)) 
	    	call error (TABLE_ACCESS, "Error loading diagnostics table")
	}

	# Fetch pointers for Object_ID and Region columns
	call clgstr ("object_col", colname, SZ_COLNAME)
	call tbcfnd (ftp, colname, fid_ptr, 1)
	if (!one_intable) 
	    call tbcfnd (dtp, colname, did_ptr, 1)

	call clgstr ("region_col", colname, SZ_COLNAME)
	call tbcfnd (ftp, colname, freg_ptr, 1)
	if (!one_intable) 
	    call tbcfnd (dtp, colname, dreg_ptr, 1)

	if (fid_ptr == NULL || did_ptr == NULL)
	    call error (BAD_OBJ_ID, "Missing Object field in input table(s)")

	# Loop through list of objects to calculate ionic abundances. 
	n_rows = tbpsta (dtp, TBL_NROWS)
	call strcpy ("Error in abundance calculation for row: %d\n", errmsg, 
			SZ_LINE)

	do drow = 1, n_rows {

	    # Determine corresponding row numbers in input/output tables. 
	    if (one_intable)
		frow = drow
	    else {
		call tbegtt (dtp, did_ptr,  drow, object, SZ_FNAME)
		call tbegtt (dtp, dreg_ptr, drow, region, SZ_FNAME)
	    	frow = match_object (ftp, fid_ptr, freg_ptr, object, region, false)
		if (frow <= 0) 
		    call error (BAD_OBJ_ID, 
				"Object_ID mis-match in input tables")
	    }

	    # Calculate dereddened line fluxes & store in structure. 
	    lfl = lfl_alloc (INDEFI)
	    iferr ( call f_ratios2 (ftp, frow, lfl) ) {
		call eprintf (errmsg)
		    call pargi (frow)
		call erract (EA_WARN)
		next
	    }

	    # Use default N_e and T_e if "constant" flag is set;
	    # otherwise, use zone-specific densities/temperatures.
	    call get_zone_nt (dtp, drow, dens, temp)

	    # Calculate ionic abundances in low-, med-, and high-ionization 
	    # zones. 
	    iferr {
		call abund2 (lfl, temp, dens) 
		call ab_out (lfl, dtp, drow)

	    } then {
		call eprintf (errmsg)
		    call pargi (frow)
		call erract (EA_WARN)
		next
	    }
	    call lfl_free (lfl)
	}

	# Close I/O tables and free memory
	if (!one_intable) 
	    call tbtclo (dtp)
	call tbtclo (ftp)
end


#-------------------------------------------------------------------------------
#  FIND_TN -	Fetch column pointers for T_e & N_e diagnostics for each 
#		zone from input table.   

int procedure find_tn (dtp, dcptr, tcptr)

#  Calling arguments
pointer	dtp			# diag table descriptor
pointer	dcptr[N_ZONES]		# density column pointers from input table
pointer	tcptr[N_ZONES]		# temperature column pointers from table

#  Local variables
char	colname[SZ_COLNAME,N_ZONES]	# names of columns
int	i			# generic
int	status			# return status: 0 == OK

begin
	status = 0

	# Fetch the column pointers for density. 
	call clgstr ("ne_low_col", colname[1,1], SZ_COLNAME)
	call clgstr ("ne_med_col", colname[1,2], SZ_COLNAME)
	call clgstr ("ne_hi_col",  colname[1,3], SZ_COLNAME)
	call tbcfnd (dtp, colname, dcptr, N_ZONES)

	# Ensure that columns exist. 
	do i = 1, N_ZONES {
	    if (dcptr[i] == NULL) 
	    	status = 1
	}

	# Fetch the column pointers for temperature. 
	call clgstr ("te_low_col", colname[1,1], SZ_COLNAME)
	call clgstr ("te_med_col", colname[1,2], SZ_COLNAME)
	call clgstr ("te_hi_col",  colname[1,3], SZ_COLNAME)
	call tbcfnd (dtp, colname, tcptr, N_ZONES)
	do i = 1, N_ZONES {
	    if (tcptr[i] == NULL) 
	    	status = 1
	}

	return (status)
end


#-------------------------------------------------------------------------------
#  GET_ZONE_NT - Fetch T_e & N_e diagnostics for each zone from input table. 

procedure get_zone_nt (tp, row, dens, temp)

#  Calling arguments
pointer	tp			# I: input table descriptor
int	row			# I: table row
real	dens[N_ZONES]		# O: density in low, med, & high ionization zone
real	temp[N_ZONES]		# O: temperature in low, med, & high ionization zone

#  Local variables
bool	const			# use constant N_e & T_e for all zones?
pointer	dcptr[N_ZONES]		# density column pointers from input table
real	N_e, T_e		# electron density, temperature from CL
int	i			# generic
bool	nulflag[N_ZONES]	# flag for INDEF values
int	status			# return status of find_tn routine
pointer	tcptr[N_ZONES]		# temperature column pointers from table

#  Functions used:
bool	clgetb()		# fetch boolean parameter from CL
real	clgetr()		# fetch real parameter from CL
int	find_tn()		# fetch T_e & N_e from diagnostics table

begin
	# Use default N_e and T_e if "constant" flag is set; 
	# otherwise, fetch diag column pointers. 
	const = clgetb ("constant")
	if (const) {
	    N_e = clgetr ("n_e")
	    T_e = clgetr ("t_e")
	    do i = 1, N_ZONES {
	    	dens[i] = N_e
	    	temp[i] = T_e
	    }
	    return

	} else {
	    status = find_tn (tp, dcptr, tcptr) 
	    if (status != OK)
		call error (2, "No Te_Low or Ne_Low in diagnostic table")
	}

	do i = 1, N_ZONES {
	    dens[i] = INDEFR
	    temp[i] = INDEFR
	}

	# Fetch T_e & N_e from columns; ensure valid values 
	# exist for low ionization zone.
	call tbrgtr (tp, dcptr, dens, nulflag, N_ZONES, row)
	call tbrgtr (tp, tcptr, temp, nulflag, N_ZONES, row)

	# Handle case of missing diagnostics for medium I.P. zones.
	if (IS_INDEF(temp[2]))
	    temp[2] = temp[1]
	if (IS_INDEF(dens[2]))
	    dens[2] = dens[1]

	# Handle case of missing diagnostics for high I.P. zones.
	if (IS_INDEF(temp[3]))
	    temp[3] = temp[2]
	if (IS_INDEF(dens[3]))
	    dens[3] = dens[2]

end



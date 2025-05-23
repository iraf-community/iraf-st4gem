# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<tbset.h>
include	<error.h>
include	"../fivel.h"
include	"../neberr.h"
include	"../neb_io.h"
include	"../zones.h"

define	DEBUG	false

#--------------------------------------------------------------------06 Feb 98--
.help t_zones.x Mar96 nebular/fivel
.ih
NAME
.nf
   t_zones - SPP driver: performs table I/O and fetches task parameters
.fi
.endhelp
#-------------------------------------------------------------------------------
#  T_ZONES -	Calculate electron temperature & densities within each 
#		of 3 ionization zones.  Performs table I/O and fetches 
#		task parameters.  

procedure t_zones ()

#  Local variables
pointer	colptr			# column descriptors from input table
pointer	otcptr[N_OTCOLS]	# column descriptors for output table
char	fluxtab[SZ_FNAME] 	# name for table of line fluxes
int	ic			# word index
pointer	ilr			# pointer to line ratio structure
pointer	itp, otp		# input, output table descriptors
char	obj_name[SZ_FNAME]	# object name to match
int	nrows			# no. table rows
int	nobj			# no. objects in input list
char	objects[SZ_FNAME]	# list of object names
char	outtab[SZ_FNAME] 	# output table of physical parameters
char	region[SZ_FNAME]	# nebula region name
pointer	reg_ptr			# pointer to region column
int	row			# input table row 
pointer	zn			# zones data structure

#  Functions used:
pointer	atl_alloc() 		# allocate AT list 
int	mat_object()		# 
bool	streq()			# are two strings equal?
int	tbpsta()		# get table header parameter
pointer	tbtopn()		# initialize table structure
int	word_fetch()		# get word from string (here, a CL parameter) 

errchk	l_ratios, obj_idcol, mat_object, zone1, zone2, zone3, z_init, z_out

begin
	# Get the input/output table names. 
	call clgstr ("fluxtab",   fluxtab, SZ_FNAME)
	call imgcluster (fluxtab, fluxtab, SZ_FNAME)
	call clgstr ("outtab",    outtab, SZ_FNAME)
	call imgcluster (outtab,  outtab, SZ_FNAME)

	# Open the input table of line fluxes; get table dimensions. 
	iferr (itp = tbtopn (fluxtab, READ_WRITE, NULL)) 
	    call error (TABLE_ACCESS, "Error loading input table")
	nrows = tbpsta (itp, TBL_NROWS)

	# Check to see if input table = output table. 
	if (streq (fluxtab, outtab)) {
	    otp = itp
	    call obj_idcol (otp, otcptr)

	# Otherwise, open/create the output table. 
	} else {
	    iferr (otp = tbtopn (outtab, NEW_FILE, NULL)) 
	    	call error (TABLE_ACCESS, "Error initializing output table")
	    call obj_idcol (otp, otcptr)
	    iferr (call tbtcre (otp))
	    	call error (TABLE_ACCESS, "Error creating output table")
	}

	call tbcfnd (itp, "Object_ID", colptr, 1)
	if (colptr == NULL) 
	    call error (BAD_OBJ_ID, "Object_ID column missing in input table")

	call tbcfnd (itp, "Region", reg_ptr, 1)

	# Create data structures for input line ratios, and for 
	# temperature/density diagnostics. 
	call malloc (ilr, LEN_ILR,   TY_STRUCT)
	call malloc (zn,  LEN_ZONES, TY_STRUCT)

	# Loop through list of objects, calculating physical parameters for 
	# each.
	call clgstr ("objects", objects, SZ_FNAME)
	nobj = 0
	ic   = 1
	while (word_fetch (objects, ic, obj_name, SZ_FNAME) > 0) {

	    # Find match between selected object and input table. 
	    row = mat_object (itp, colptr, reg_ptr, obj_name, region, true) 
	    if (row > 0) {
		nobj = nobj + 1
		call tbeptt (otp, otcptr[ID_COL], nobj, obj_name)
		call tbeptt (otp, otcptr[RG_COL], nobj, region)

		# Calculate nebular physical parameters in each of three zones.
		iferr {
		    Z_ATL(zn) = atl_alloc (INDEFI)

		    # Calculate line ratios. 
		    call l_ratios (itp, ilr, row)

		    call z_init (zn)

		    call zone1 (zn, ilr)
		    call zone2 (zn, ilr)
		    call zone3 (zn, ilr)
		    call z_out (zn, otp, nobj)
		    call atl_free (Z_ATL(zn), true)

		} then {
		    call erract (EA_WARN)
		    call eprintf ("Error in T_e/N_e calculations for object: %s, region: %s\n")
			call pargstr (obj_name)
			call pargstr (region)
		    call flush (STDERR)
		    #next
		}
	    }
	}

	#if (nobj == 0) 
	    #call error (BAD_OBJ_ID, "No matching Object_ID in input table")

	# Close I/O tables and free memory. 
	if (itp != otp) 
	    call tbtclo (otp)
	call tbtclo (itp)
	call mfree (ilr, TY_STRUCT)
	call atl_free (Z_ATL(zn), true)
	call mfree (zn, TY_STRUCT)
end



# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<tbset.h>
include	<error.h>
include	"../fivel.h"
include	"../neb_io.h"
include	"../zones.h"

define	DEBUG	false

#---------------------------------------------------------------------5 Feb 97--
.help z_out.x Mar96 nebular/fivel
.ih
NAME
.nf
    z_out - Write N_e, T_e diagnostics to output table
obj_idcol - Allocate the output table columns for object/region ID.  
  diagcol - Allocate the output table columns  
.fi
.endhelp
#-------------------------------------------------------------------------------
#  Z_OUT --	Write N_e, T_e diagnostics to output table.  

procedure z_out (zn, tp, row)

#  Calling arguments:
pointer	zn			# I: zones data structure
pointer	tp			# I: table descriptor
int	row			# I: output table row number

#  Local variables:
pointer	cp[N_OTCOLS]		# I: table column descriptors

begin
	call diagcol (tp, cp)

	# Write low ionization diagnostics to output table.
	call tbeptr (tp, cp[N_O2_COL], row, Ne_Oii(zn))
	call tbeptr (tp, cp[N_S2_COL], row, Ne_Sii(zn))

	call tbeptr (tp, cp[T_N2_COL], row, Te_Nii(zn))
	call tbeptr (tp, cp[T_O2_COL], row, Te_Oii(zn))
	call tbeptr (tp, cp[T_S2_COL], row, Te_Sii(zn))

	# Write medium ionization diagnostics to output table.
	call tbeptr (tp, cp[N_CL3_COL], row, Ne_CLiii(zn))
	call tbeptr (tp, cp[N_AR4_COL], row, Ne_ARiv(zn))
	call tbeptr (tp, cp[N_C3_COL],  row, Ne_Ciii(zn))

	call tbeptr (tp, cp[T_O3_COL],  row, Te_Oiii(zn))
	call tbeptr (tp, cp[T_NE3_COL], row, Te_NEiii(zn))
	call tbeptr (tp, cp[T_AR3_COL], row, Te_ARiii(zn))
	call tbeptr (tp, cp[T_AR4_COL], row, Te_ARiv(zn))
	call tbeptr (tp, cp[T_S3_COL],  row, Te_Siii(zn))

	# Write high ionization diagnostics to output table. 
	call tbeptr (tp, cp[N_NE4_COL], row, Ne_NEiv(zn))

#	call tbeptr (tp, cp[T_NE4_COL], row, Te_NEiv(zn))
	call tbeptr (tp, cp[T_NE5_COL], row, Te_NEv(zn))
	call tbeptr (tp, cp[T_AR5_COL], row, Te_ARv(zn))

	# Write summary diagnostics to output table. 
	call tbeptr (tp, cp[N_LOW_COL], row, Ne_Low(zn))
	call tbeptr (tp, cp[T_LOW_COL], row, Te_Low(zn))
	call tbeptr (tp, cp[N_MED_COL], row, Ne_Med(zn))
	call tbeptr (tp, cp[T_MED_COL], row, Te_Med(zn))
	call tbeptr (tp, cp[N_HI_COL],  row, Ne_Hi(zn))
	call tbeptr (tp, cp[T_HI_COL],  row, Te_Hi(zn))
end


#-------------------------------------------------------------------------------
#  OBJ_IDCOL --	Allocate the output table columns for object/region ID.  

procedure obj_idcol (tp, otcp)

#  Calling arguments:
pointer	tp			# Table descriptor
pointer	otcp[ARB]		# Output column descriptors

#  Local variables:
int	col			# loop counter
char	colname[SZ_COLNAME+1,2]	# Column names
char	colunits[SZ_COLUNITS+1,2] # Column units
char	colfmt[SZ_COLFMT+1,2]	# Column formats
int	cdtype[2]		# Column data types
int	lendata			# Character sizes

begin
	call strcpy ("Object_ID", colname[1,ID_COL],  SZ_COLNAME)
	call strcpy ("none",      colunits[1,ID_COL], SZ_COLUNITS)
	call strcpy ("%14s",      colfmt[1,ID_COL],   SZ_COLFMT)
	cdtype[ID_COL] = -14

	call strcpy ("Region", colname[1,RG_COL],  SZ_COLNAME)
	call strcpy ("none",   colunits[1,RG_COL], SZ_COLUNITS)
	call strcpy ("%2s",    colfmt[1,RG_COL],   SZ_COLFMT)
	cdtype[RG_COL] = -2

	# Define these table columns if non-existant
	call tbcfnd (tp, colname, otcp, 2)
	do col = 1, 2 {
	    if (otcp[col] == NULL) 
		call tbcdef (tp, otcp[col], colname[1,col], colunits[1,col],
			colfmt[1,col], cdtype[col], lendata, 1)
	}
end


#-------------------------------------------------------------------------------
#  DIAGCOL --	Allocate the output table columns.  If a column already 
#		exists, use it, assuming it's attributes are appropriate.  

procedure diagcol (tp, otcp)

#  Calling arguments:
pointer	tp				# Table descriptor
pointer	otcp[N_OTCOLS]			# Output column descriptors

#  Local variables:
int	col				# loop counter
char	colname[SZ_COLNAME+1,N_OTCOLS]	# Column names
char	colunits[SZ_COLUNITS+1,N_OTCOLS] # Column units
char	colfmt[SZ_COLFMT+1,N_OTCOLS]	# Column formats
int	cdtype[N_OTCOLS]		# Column data types
int	lendata				# Character sizes

begin
	# Output column names. 
	call strcpy ("Ne_Ciii]",   colname[1,N_C3_COL],  SZ_COLNAME)
	call strcpy ("Ne_[Oii]",   colname[1,N_O2_COL],  SZ_COLNAME)
	call strcpy ("Ne_[NEiv]",  colname[1,N_NE4_COL], SZ_COLNAME)
	call strcpy ("Ne_[Sii]",   colname[1,N_S2_COL],  SZ_COLNAME)
	call strcpy ("Ne_[CLiii]", colname[1,N_CL3_COL], SZ_COLNAME)
	call strcpy ("Ne_[ARiv]",  colname[1,N_AR4_COL], SZ_COLNAME)
	call strcpy ("Te_[Nii]",   colname[1,T_N2_COL],  SZ_COLNAME)
	call strcpy ("Te_[Oii]",   colname[1,T_O2_COL],  SZ_COLNAME)
	call strcpy ("Te_[Oiii]",  colname[1,T_O3_COL],  SZ_COLNAME)
	call strcpy ("Te_[NEiii]", colname[1,T_NE3_COL], SZ_COLNAME)
#	call strcpy ("Te_[NEiv]",  colname[1,T_NE4_COL], SZ_COLNAME)
	call strcpy ("Te_[NEv]",   colname[1,T_NE5_COL], SZ_COLNAME)
	call strcpy ("Te_[Sii]",   colname[1,T_S2_COL],  SZ_COLNAME)
	call strcpy ("Te_[Siii]",  colname[1,T_S3_COL],  SZ_COLNAME)
	call strcpy ("Te_[ARiii]", colname[1,T_AR3_COL], SZ_COLNAME)
	call strcpy ("Te_[ARiv]",  colname[1,T_AR4_COL], SZ_COLNAME)
	call strcpy ("Te_[ARv]",   colname[1,T_AR5_COL], SZ_COLNAME)
	call strcpy ("Ne_Low",     colname[1,N_LOW_COL], SZ_COLNAME)
	call strcpy ("Ne_Med",     colname[1,N_MED_COL], SZ_COLNAME)
	call strcpy ("Ne_Hi",      colname[1,N_HI_COL],  SZ_COLNAME)
	call strcpy ("Te_Low",     colname[1,T_LOW_COL], SZ_COLNAME)
	call strcpy ("Te_Med",     colname[1,T_MED_COL], SZ_COLNAME)
	call strcpy ("Te_Hi",      colname[1,T_HI_COL],  SZ_COLNAME)

	# Output column units
	call strcpy ("N/cm^3", colunits[1,N_C3_COL],  SZ_COLUNITS)
	call strcpy ("N/cm^3", colunits[1,N_O2_COL],  SZ_COLUNITS)
	call strcpy ("N/cm^3", colunits[1,N_NE4_COL], SZ_COLUNITS)
	call strcpy ("N/cm^3", colunits[1,N_S2_COL],  SZ_COLUNITS)
	call strcpy ("N/cm^3", colunits[1,N_CL3_COL], SZ_COLUNITS)
	call strcpy ("N/cm^3", colunits[1,N_AR4_COL], SZ_COLUNITS)
	call strcpy ("Kelvin", colunits[1,T_N2_COL],  SZ_COLUNITS)
	call strcpy ("Kelvin", colunits[1,T_O2_COL],  SZ_COLUNITS)
	call strcpy ("Kelvin", colunits[1,T_O3_COL],  SZ_COLUNITS)
	call strcpy ("Kelvin", colunits[1,T_NE3_COL], SZ_COLUNITS)
#	call strcpy ("Kelvin", colunits[1,T_NE4_COL], SZ_COLUNITS)
	call strcpy ("Kelvin", colunits[1,T_NE5_COL], SZ_COLUNITS)
	call strcpy ("Kelvin", colunits[1,T_S2_COL],  SZ_COLUNITS)
	call strcpy ("Kelvin", colunits[1,T_S3_COL],  SZ_COLUNITS)
	call strcpy ("Kelvin", colunits[1,T_AR3_COL], SZ_COLUNITS)
	call strcpy ("Kelvin", colunits[1,T_AR4_COL], SZ_COLUNITS)
	call strcpy ("Kelvin", colunits[1,T_AR5_COL], SZ_COLUNITS)

	call strcpy ("N/cm^3", colunits[1,N_LOW_COL], SZ_COLUNITS)
	call strcpy ("N/cm^3", colunits[1,N_MED_COL], SZ_COLUNITS)
	call strcpy ("N/cm^3", colunits[1,N_HI_COL],  SZ_COLUNITS)
	call strcpy ("Kelvin", colunits[1,T_LOW_COL], SZ_COLUNITS)
	call strcpy ("Kelvin", colunits[1,T_MED_COL], SZ_COLUNITS)
	call strcpy ("Kelvin", colunits[1,T_HI_COL],  SZ_COLUNITS)

	# Set the column output formats, data types, & char column sizes. 
	lendata = 0

	# Loop over the remaining columns. 
	do col = 3, N_OTCOLS {
	    call strcpy ("%8.0f", colfmt[1,col], SZ_COLFMT)
	    cdtype[col] = TY_REAL
	}

	# Check for existing columns.  We won't worry for now about columns 
	# that don't match in datatype. 

	call tbcfnd (tp, colname, otcp, N_OTCOLS)

	do col = 3, N_OTCOLS {
	    if (otcp[col] == NULL)	#  Define output table columns 

		call tbcdef (tp, otcp[col], colname[1,col], colunits[1,col],
			colfmt[1,col], cdtype[col], lendata, 1)
	}
end



# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<tbset.h>
include	"../flux.h"

#---------------------------------------------------------------------9 Jul 97--
.help ab_col.x Jul97 nebular/abund
.ih
NAME
.nf
   ab_out - Output ionic abundances to table row
  ab_cols - Return all column attributes for the ionic abundance 
.fi
.endhelp
#-------------------------------------------------------------------------------
# AB_OUT -	Output ionic abundances to table row. 

procedure ab_out (lfl, otp, row)

#  Calling arguments:
pointer	lfl			# I: flux list data structure
pointer	otp			# I: output table descriptor
int	row			# I: output table row

#  Declarations:
int	cdtype			# Column data type
char	format[SZ_COLFMT]	# Column format
pointer	lf			# flux data structure
int	i			# generic
int	n_ions			# number of ions for which to compute abundance
char	name[SZ_COLNAME]	# Column name
pointer	otcp			# output table column descriptors 
char	units[SZ_COLUNITS]	# Column units

errchk	ab_cols

begin
	n_ions = LFL_N(lfl)
	do i = 1, n_ions {
	    lf = LFL_A(lfl, i)
	    call ab_cols (LF_ATOM(lf), LF_ION(lf), name, units, format, cdtype)

	    # Check for existing column name and either create a new column...
	    call tbcfnd (otp, name, otcp, 1)
	    if (otcp == NULL)
		call tbcdef (otp, otcp, name, units, format, cdtype, 0, 1)

	    # ...or use the existing one.
	    else {

		# Ensure that existing column has correct units & format. 
	    	call tbcnit (otp, otcp, units)
	    	call tbcfmt (otp, otcp, format)
	    }

	    # Write the ionic abundance
	    call tbeptr (otp, otcp, row, LF_ABUND(lf))
	}
end


#-------------------------------------------------------------------------------
# AB_COLS -	Return all column attributes for the ionic abundance. 

procedure ab_cols (atom, ion, name, units, format, cdtype)

#  Calling arguments:
int	atom			# I: atomic number
int	ion			# I: ionization stage
char	name[SZ_COLNAME]	# O: Column name
char	units[SZ_COLUNITS]	# O: Column units
char	format[SZ_COLFMT]	# O: Column format
int	cdtype			# O: Column data type

#  Declarations:
string	s1 	"Ni_(%s^+%d)"
char	symbol[3]		# atomic symbol

errchk	get_atomic_symbol

begin
	# Construct the column name
	call get_atomic_symbol (atom, symbol, 3)
	call sprintf (name, SZ_COLNAME, s1)
	    call pargstr (symbol)
	    call pargi (ion)

	call strcpy ("1/cm^3", units, SZ_COLUNITS)
	call strcpy ("%8.4g", format, SZ_COLFMT)
	cdtype = TY_REAL
end



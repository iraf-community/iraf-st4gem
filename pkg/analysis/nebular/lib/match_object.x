# Copyright restrictions apply - see st4gem$copyright.stsdas

include	<tbset.h>
include	"../neberr.h"

define	WILDCARD	"*"
define	NO_MATCH	0

#---------------------------------------------------------------------7 Sep 97--
.help match_object.x Feb97 nebular/fivel
.ih
NAME
.nf
match_object - Match the object & region ID strings to input table entries
   obj_match - Determine whether a match on the object/region occured
  mat_object - Match the object and region ID strings to input table entries
.fi
.ih
DESCRIPTION
Match the object and region ID strings to input table entries.  
Returns the table row number of the match, or zero if no match 
is found.
.endhelp
#-------------------------------------------------------------------------------
#  MATCH_OBJECT -	Match the object & region ID strings to input table 
#			entries.  

int procedure match_object (itp, obj_col, reg_col, object, region, wild_ok)

#  Calling arguments:
pointer	itp			# I: input table descriptor
pointer	obj_col			# I: Object column pointer
pointer	reg_col			# I: Region column pointer
char	object[ARB]		# I: object name to match
char	region[ARB]		# I: region name to match
bool	wild_ok			# I: wildcard match permitted?
int	row			# O: table row number of match

#  Declarations:
int	n_rows			# number of table rows
bool	obj_match()		# did the target match the object/region ID?
char	reg_ID[SZ_FNAME]	# matched region name from input table
char	targ_ID[SZ_FNAME]	# matched Object_ID name from input table
int	tbpsta()		# get table header parameter

begin
	n_rows = tbpsta (itp, TBL_NROWS)
	if (obj_col == NULL)
	    return (NO_MATCH)

	for (row = 1; row <= n_rows; row = row+1) {
	    call tbegtt (itp, obj_col, row, targ_ID, SZ_FNAME)
	    if ( obj_match (object, targ_ID, SZ_FNAME, wild_ok ) ) {
		if (reg_col == NULL) 
		    return (row)
		else {
		    call tbegtt (itp, reg_col, row, reg_ID, SZ_FNAME)
		    if ( obj_match (region, reg_ID, SZ_FNAME, wild_ok ) ) 
			return (row)
		}
	    }
	}

	return (NO_MATCH)
end


#-------------------------------------------------------------------------------
#  OBJ_MATCH -	Determine whether a match on the object/region occured. 
#		If the match is to a wildcard, replace with the target name. 

bool procedure obj_match (object, target, sz_name, wild_ok)

#  Calling arguments:
char	object[sz_name]		# I/O: object name to match
char	target[sz_name]		# I: matched Object_ID name from input table
int	sz_name			# I: size of object/target names
bool	wild_ok			# I: wildcard match permitted?

#  Declarations:
bool	streq()			# are two strings equal?

begin
	if ( wild_ok && streq (object[1], WILDCARD) ) {
	    call strcpy (target, object, sz_name)
	    return (true)

	} else if ( streq (object, target) )
	    return (true)

	else
	    return (false)
end


#-------------------------------------------------------------------------------
#  MAT_OBJECT -	Match the object and region ID strings to input table 
#			entries.  

int procedure mat_object (itp, obj_col, reg_col, object, region, wild_ok)

#  Calling arguments:
pointer	itp			# I: input table descriptor
pointer	obj_col			# I: Object column pointer
pointer	reg_col			# I: Region column pointer
char	object[ARB]		# I: object name to match
char	region[ARB]		# I: matched region name from input table
bool	wild_ok			# I: wildcard match permitted?
int	row			# O: table row number of match

#  Declarations:
int	n_rows			# I: number of table rows
bool	obj_match()		# did the target match the object/region ID?
#char	reg_ID[SZ_FNAME]	# matched region name from input table
char	targ_ID[SZ_FNAME]	# matched Object_ID name from input table
#bool	streq()			# are two strings equal?
int	tbpsta()		# get table header parameter

begin
	n_rows = tbpsta (itp, TBL_NROWS)
	if (obj_col == NULL)
	    return (NO_MATCH)

	for (row = 1; row <= n_rows; row = row+1) {
	    call tbegtt (itp, obj_col, row, targ_ID, SZ_FNAME)
	    if ( obj_match (object, targ_ID, SZ_FNAME, wild_ok) ) {
		if (reg_col == NULL) 
		    return (row)
		else {
		    call tbegtt (itp, reg_col, row, region, SZ_FNAME)
		    return (row)
		}
	    }
	}

	return (NO_MATCH)
end



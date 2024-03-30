# Copyright restrictions apply - see st4gem$copyright.stsdas 

include	<tbset.h>
include <synphot.h>
include "libsynphot.h"

#* HISTORY *
#* B.Simon	20-Oct-97	Check for parameterized columns in tabrange
#* B.Simon	02-Feb-98	Fixed bug in z function
#* B.Simon	18-Sep-00	Rewritten to call synwave

# CALCRANGE -- Calculate wavelength range appropriate for expression

procedure calcrange (pcode, maxcode, graphtab, comptab, minwave, maxwave)

int	pcode[ARB]	# i: pseudocode used by calculator
int	maxcode		# i: maximum length of pseudocode
char    graphtab[ARB]   # i: graph table name
char    comptab[ARB]    # i: component lookup table name
real	minwave		# o: short end of wavelength range
real	maxwave		# o: long end of wavelength range
#--
int	status
pointer	rg1,rg2,rg3

string	badvega  "Error in Vega spectrum"

pointer	syw_calc(), syw_create(), syw_intersect()

begin
	# Calculate the range set

	rg1 = syw_calc (pcode, maxcode, graphtab, comptab)
	call syw_bound (rg1, DEFMIN, DEFMAX)

	# Replace INDEF values with defaults

	call syw_bound (rg1, DEFMIN, DEFMAX)

	# Take intersection with vega's spectrum

	call tabrange (VEGA, minwave, maxwave, status)
	if (status == ERR)
	    call synphoterr (badvega, VEGA)

	rg2 = syw_create (minwave, maxwave)
	rg3 = syw_intersect (rg1, rg2)

	call syw_free (rg1)
	call syw_free (rg2)

	# Find minimum and maximum range values

	call syw_limits (rg3, minwave, maxwave)
	call syw_free (rg3)

end

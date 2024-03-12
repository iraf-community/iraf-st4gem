###########################################################################
#                    Center for Astrophysical Sciences
#                        Johns Hopkins University
#
#  Synopsis:	sfinit()
#
#  Description:	SFINIT intializes data in global common for SPECFIT
#
#  Arguments:	none
#
#  Returns:	none
#
#  Notes:	Information shared in common blocks defined in "specfit.com".
#
#  History:	May 1989	Gerard Kriss
#		August 1989	gak	changed empirical to 2 parameters
#		March 1989	gak	Added recomb and multigauss
#		Sep 30 1991	gak	Added initialization for usercont
#		Oct  1 1991	gak	Same for userline and userabs
#		Oct  1 1991	gak	made default "nofluxes" for /ansr/
#		August, 94	J Grimes Created sfinitcomp for the new tasks
#
###########################################################################

include	"specfit.h"

procedure sfinit()

int i

include	"specfit.com"

begin

	call sfinitcomp()	

	cfirst = 1
	pfirst = 1
	afirst = 1

	ncont = -1	# If no flux interval database file is specified,
			# the call to sfgetflux will be skipped.

#	This call will invoke a routine to ignore floating point errors, but
#	it will only compile under SunOS4.x
#
#	call igfpe()

	for ( i = 1; i<=NCOMP; i = i + 1 ) {
		inplot[i] = YES
		}

end




























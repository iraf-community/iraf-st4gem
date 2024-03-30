include "libsynphot.h"

# GRAFFILES -- Return an array of component throughput table names
#
# This routine merely returns the results of searchgraf minus any
# parameters that mave have been associated with it.  It is used to 
# provide information about the workings of synphot.

procedure graffiles (mode, graphtab, comptab, maxname, maxthru, nthru, thruput)

char	mode[ARB]			# i: instrument mode
char	graphtab[ARB]			# i: graph table name
char	comptab[ARB]			# i: component lookup table name
int	maxname				# i: max length of thruput table name
int	maxthru				# i: max number of throughput tables
int	nthru				# o: actual number of throughput tables
char	thruput[maxname,maxthru]	# o: throughput table names
#--
bool	verbose
pointer	sp, nparam, paramlist

data	verbose	 / false /

errchk	searchgraf

begin
	# Allocate arrays for parameter variables that we will
	# drop on the floor

	call smark (sp)
	call salloc (nparam, maxthru, TY_INT)
	call salloc (paramlist, MAXPARAM*maxthru, TY_REAL)

	# Search graph table for list of throughput table names
	# and their associated parameters

#	call searchgraf (verbose, graphtab, 
        call searchgraf (verbose, graphtab, GRF_COMPID,   # Calling sequence change VGL 9/28/01
                         comptab, mode, maxthru, 
			 MAXPARAM, maxname, nthru, Memi[nparam], 
			 Memr[paramlist], thruput)

	call sfree (sp)
end

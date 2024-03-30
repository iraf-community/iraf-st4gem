include "libsynphot.h"

# BANDNAME -- Return names of throughput tables used by bandpass function

procedure bandname (obsmode, graphtab, comptab, log)

char	obsmode[ARB]	# i: observation mode string
char	graphtab[ARB]	# i: graph table name
char	comptab[ARB]	# i: component lookup table name
int	log		# i: log file names are written on
#--
int	icomp, ncomp
pointer	sp, flist, fname

errchk	graffiles

begin
	# Allocate memory for temporarary arrays

	call smark (sp)
	call salloc (flist, (SZ_FNAME+1)*MAXLIST, TY_CHAR)

	# Produce list of table names from the observation mode

	call graffiles (obsmode, graphtab, comptab, SZ_FNAME, 
			MAXLIST, ncomp, Memc[flist])

	# Print the list to the log file

	fname = flist
	do icomp = 1, ncomp {
	    call fprintf (log, "%s\n")
	    call pargstr (Memc[fname])

	    fname = fname + SZ_FNAME + 1
	}

	call sfree (sp)
end

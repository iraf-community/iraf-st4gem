include "libsynphot.h"

# LISTPATH -- Create list of component throughput tables

procedure listpath (mode, graphtab, comptab, path, mxpath)

char	mode[ARB]	# i: instrument mode
char	graphtab[ARB]	# i: graph table name
char	comptab[ARB]	# i: component lookup table name
char	path[ARB]	# o: list of component throughput tables
int	mxpath		# i: maximum length of path string
#--
int	ncomp
pointer	sp, filelist

errchk	graffiles

begin
	# Allocate memory for temporary arrays

	call smark (sp)
	call salloc (filelist, (SZ_FNAME+1)*MAXLIST, TY_CHAR)

	# Get list of component throughput files

	call inisyntab

	call graffiles (mode, graphtab, comptab, SZ_FNAME, 
			MAXLIST, ncomp, Memc[filelist])

	call clssyntab

	# Concatenate file names into path string

	call wpath (SZ_FNAME, ncomp, Memc[filelist], mxpath, path)

	call sfree (sp)
end

include "libsynphot.h"

# FINDTHRUPUT -- Find thruput table containing boolean keyword = T

procedure findthruput (keyword, obsmode, graphtab, comptab, tabname, maxch)

char	keyword[ARB]	# i: keyword name
char	obsmode[ARB]	# i: observation mode string
char	graphtab[ARB]	# i: graph table name
char	comptab[ARB]	# i: component lookup table name
char	tabname[ARB]	# o: thruput table name
int	maxch		# i: maximum length of table name
#--
bool	match, value
int	icomp, ncomp
pointer	sp, tp, flist, fname

string	toomany  "findthruput: more than one table contains keyword"

bool	rdtabhdb() 
pointer	opnsyntab()

errchk	graffiles

begin
	# Allocate memory for temporarary arrays

	call smark (sp)
	call salloc (flist, (SZ_FNAME+1)*MAXLIST, TY_CHAR)

	# Produce list of table names from the observation mode

	call graffiles (obsmode, graphtab, comptab, SZ_FNAME, 
			MAXLIST, ncomp, Memc[flist])

	# Search the list of thruput table names for the keyword

	match = false
	fname = flist
	do icomp = 1, ncomp {
	    tp = opnsyntab (Memc[fname])

	    ifnoerr {
		value = rdtabhdb (tp, keyword)

	    } then {
		if (value) {
		    if (match)
			call synphoterr (toomany, Memc[fname])

		    match = true
		    call strcpy (Memc[fname], tabname, maxch)
		}
	    }

	    fname = fname + SZ_FNAME + 1
	}

	# Null string if no match

	if (! match)
	    tabname[1] = EOS

	call sfree (sp)
end

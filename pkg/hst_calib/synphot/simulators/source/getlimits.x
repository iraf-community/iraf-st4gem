#* HISTORY *
#* B.Simon	15-Aug-95	original

include "limit.h"

# GETLIMITS -- Read count rate limits from catalog

procedure getlimits (limitcat, obsmode, limit)

char	limitcat[ARB]	# i: limit catalog
char	obsmode[ARB]	# i: observation mode
real	limit[LEN_L]	# o: limit parameters
#--
int	jrow, icol, ic
pointer	sp, catalog, cname, tp, cp

string	obscol     "OBSMODE"
string	limitcols  LIMIT_STR  

pointer	tbtopn(), word_fetch()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (catalog, SZ_FNAME, TY_CHAR)
	call salloc (cname, SZ_FNAME, TY_CHAR)

	# Open aperture dimension  catalog

	call lastfile (limitcat, Memc[catalog], SZ_FNAME)
	tp = tbtopn (Memc[catalog], READ_ONLY, 0)

	# Find table row which matches obsmode

	call findmode (tp, obscol, obsmode, jrow)

	# Copy limit parameters from table

	if (jrow == 0) {
	    call amovkr (INDEFR, limit, LEN_L)

	} else {
	    ic = 1
	    icol = 1
	    while (word_fetch (limitcols, ic, Memc[cname], SZ_FNAME) > 0) {
		call syncolptr (tp, Memc[cname], icol+1, cp)
		call tbegtr (tp, cp, jrow, limit[icol])

		icol = icol + 1
	    }
	}

	call tbtclo (tp)
	call sfree (sp)
end


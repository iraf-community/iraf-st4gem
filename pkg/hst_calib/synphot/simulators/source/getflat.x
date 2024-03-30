include <tbset.h>

#* HISTORY *
#* B.Simon	30-Jun-95	original

# GETFLAT -- Read flat field from catalog

procedure getflat (flatcat, obsmode, flatfile, maxch)

char	flatcat[ARB]	# i: Catalog of flat fields
char	obsmode[ARB]	# i: Observation mode
char	flatfile[ARB]	# o: Flatfield filename
int	maxch		# i: Maximum length of flatfield file
#--
int	jrow
pointer	sp, catalog, tp, cp

string	obscol   "OBSMODE"
string	filecol  "FILENAME"

pointer	tbtopn()

begin
	# Return empty flat field if catalog not found

	if (flatcat[1] == EOS) {
	    flatfile[1] = EOS
	    return
	}

	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (catalog, SZ_FNAME, TY_CHAR)

	# Open aperture dimension  catalog

	call lastfile (flatcat, Memc[catalog], SZ_FNAME)
	tp = tbtopn (Memc[catalog], READ_ONLY, 0)

	# Find table row which matches obsmode

	call findmode (tp, obscol, obsmode, jrow)

	# Copy filename from table

	if (jrow == 0) {
	    flatfile[1] = EOS

	} else {
	    call syncolptr (tp, filecol, 2, cp)
	    call tbegtt (tp, cp, jrow, flatfile, maxch)
	}

	call tbtclo (tp)
	call sfree (sp)
end

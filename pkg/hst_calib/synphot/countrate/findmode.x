include	<tbset.h>

#* HISTORY *
#* B.Simon	20-Mar-95	original

# FINDMODE -- Find the best match  to an observation mode in a catalog

bool procedure findmode (catalog, obsmode, output, maxch)

char	catalog[ARB]	# i: catalog name
char	obsmode[ARB]	# i: observation mode
char	output[ARB]	# o: output string
int	maxch		# i: length of output string
#--
int	maxkey, nkey, nrow, irow, jrow, ic
pointer	sp, file, catmode, key, tp, cp

string	obscol   "OBSMODE"
string	filecol  "FILENAME"

int	tbpsta(), word_fetch(), word_match()
pointer	tbtopn()

begin
	# Allocate memory for strings

	call smark (sp)
	call salloc (file, SZ_FNAME, TY_CHAR)
	call salloc (catmode, SZ_FNAME, TY_CHAR)
	call salloc (key, SZ_FNAME, TY_CHAR)

	# Open catalog

	call lastfile (catalog, Memc[file], SZ_FNAME)
	tp = tbtopn (Memc[file], READ_ONLY, 0)
	nrow = tbpsta (tp, TBL_NROWS)

	# Search for best match with observation mode

	call syncolptr (tp, obscol, 1, cp)

	jrow = 0
	maxkey = 0
	do irow = 1, nrow {
	    call tbegtt (tp, cp, irow, Memc[catmode], SZ_FNAME)

	    # Check each keyword in catalog mode against input obsmode

	    ic = 1
	    nkey = 0
	    while (word_fetch (Memc[catmode], ic, Memc[key], SZ_FNAME) > 0) {
		if (word_match (Memc[key], obsmode) == 0) {
		    nkey = 0
		    break
		}

		nkey = nkey + 1
	    }

	    # Save best match so far

	    if (nkey > maxkey) {
		maxkey = nkey
		jrow = irow
	    }
	}

	# Copy filename associated with obsmode to output string

	if (jrow == 0) {
	    output[1] = EOS

	} else {
	    call syncolptr (tp, filecol, 2, cp)
	    call tbegtt (tp, cp, jrow, output, maxch)
	}

	# Return flag which indicates success of match

	call sfree (sp)
	return (jrow != 0)

end

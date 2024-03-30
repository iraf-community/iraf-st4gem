include	<tbset.h>

#* HISTORY *
#* B.Simon	11-Jul-95	extracted from countrate/modefile

# FINDMODE -- Find the row in a catalog table matching an obsmode

procedure findmode (tp, obscol, obsmode, jrow)

pointer	tp		# i: table descriptor
char	obscol[ARB]	# i: obsmode column name
char	obsmode[ARB]	# i: observation mode
int	jrow		# o: best matching mode
#--
int	maxkey, nkey, nrow, irow, ic
pointer	sp, catmode, key, cp

int	tbpsta(), word_fetch(), word_match()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (catmode, SZ_FNAME, TY_CHAR)
	call salloc (key, SZ_FNAME, TY_CHAR)

	# Search for best match

	jrow = 0
	maxkey = 0
	nrow = tbpsta (tp, TBL_NROWS)
	call syncolptr (tp, obscol, 1, cp)

	do irow = 1, nrow {
	    call tbegtt (tp, cp, irow, Memc[catmode], SZ_FNAME)

	    # Check each keyword in catalog obsmode against input obsmode

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

	call sfree (sp)
end


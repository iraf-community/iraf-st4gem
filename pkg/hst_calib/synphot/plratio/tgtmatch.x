include	<tbset.h>

#* HISTORY *
#* B.Simon	14-Jul-94	original

# TGTMATCH -- Find obsmode and target in photometry file that matches spectrum

procedure tgtmatch (obsmode, spectrum, vzero, pfile, tgtlist, ntgt)

char	obsmode[ARB]	# i: instrument observation mode
char	spectrum[ARB]	# i: target spectrum
char	vzero[ARB]	# i: variable list
char	pfile[ARB]	# i: photometry file
pointer	tgtlist		# o: target list
int	ntgt		# o: number of targets
#--
bool	found
int	nomode, offset, nrow, itgt, irow
pointer	sp, phottab, mode, spec, fspec, pmode, pspec, errmsg, tablist
pointer	tp, pptr, optr, sptr, obs, tgt
real	v0

string	matchfmt  "%s[%d]"
string	nomatch   "Photometry table row unmatched"

bool	strne(), streq()
int	numlist(), numvzero(), nxtlist(), nxtvzero(), tbpsta()
pointer	rdlist(), tbtopn()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (phottab, SZ_LINE, TY_CHAR)
	call salloc (mode, SZ_LINE, TY_CHAR)
	call salloc (spec, SZ_LINE, TY_CHAR)
	call salloc (fspec, SZ_LINE, TY_CHAR)
	call salloc (pmode, SZ_LINE, TY_CHAR)
	call salloc (pspec, SZ_LINE, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	# Replace empty obsmode with 1.0

	if (obsmode[1] == EOS) {
	    nomode = YES
	    call strcpy ("1.0", Memc[mode], SZ_LINE)
	} else {
	    nomode = NO
	    call strcpy (obsmode, Memc[mode], SZ_LINE)
	}

	# Initialize table cache and synphot variables

	call inisyntab
	call undefsynvar

	# Process the pfile, obsmode, spectrum, and vzero strings

	pptr = rdlist (pfile)
	optr = rdlist (Memc[mode])
	sptr = rdlist (spectrum)
	call rdvzero (vzero)

	tgtlist = NULL
	ntgt = 0

	# Read each photomery table, looking for a combination of obsmode 
	# and target that match the obsmode and spectrum input by the user
	# If a match is found, store the row number in the table plus
	# the length of any previous tables in the tgtlist array

	offset = 0
	while (nxtlist (pptr, Memc[phottab], SZ_LINE) != EOF) {
	    iferr {
		tp = tbtopn (Memc[phottab], READ_ONLY, NULL)
	    } then {
		next
	    }

	    nrow = tbpsta (tp, TBL_NROWS)
	    call syncolptr (tp, "OBSMODE", 3, obs)
	    call syncolptr (tp, "TARGETID", 4, tgt)

	    itgt = 0
	    while (nxtlist (optr, Memc[mode], SZ_LINE) != EOF) {
		while (nxtlist (sptr, Memc[spec], SZ_LINE) != EOF) {
		    while (nxtvzero (v0) != EOF) {
			# Allocate array on first pass thru loop

			if (tgtlist == NULL)  {
			    ntgt = numlist(optr) * numlist (sptr) * numvzero ()
			    call calloc (tgtlist, ntgt, TY_INT)

			    call salloc (tablist, ntgt, TY_INT)
			    call amovkr (NO, Memi[tablist], ntgt)
			}

			# Skip this combination of obsmode and spectrum if
			# a match has already been found

			if (Memi[tgtlist+itgt] != 0)
			    next

			call fillexpr (Memc[spec], Memc[fspec], SZ_LINE)

			# Loop over each row of table, looking for match

			do irow = 1, nrow {
			    call tbegtt (tp, obs, irow, Memc[pmode], SZ_LINE)
			    call tbegtt (tp, tgt, irow, Memc[pspec], SZ_LINE)

			    # Does the spectrum match?

			    if (strne (Memc[fspec], Memc[pspec]))
				next

			    # An empty obsmode string will match an obsmode
			    # that is blank, "none", or "1.0" in the table
			    
			    if (nomode == NO) {
				found = streq (Memc[mode], Memc[pmode])

			    } else {
				found = false
				if (streq (Memc[pmode], Memc[mode])) {
				    found = true
				    
				} else {
				    call strfix (Memc[pmode])
				    if (Memc[pmode] == EOS) {
					found = true
				    } else if (streq (Memc[pmode], "none")) {
					found = true
				    } else {
					found = false
				    }
				}
			    }

			    # Update target array if match was found and
			    # skip checking remaining rows of table

			    if (found) {
				Memi[tgtlist+itgt] = irow + offset
				Memi[tablist+offset+irow-1] = YES
				break
			    }
			}

			itgt = itgt + 1
		    }
		}
	    }

	    # Look for unmatched table rows and issue warning message

	    do irow = 1, nrow {
		if (Memi[tablist+offset+irow-1] == NO) {
		    call sprintf (Memc[errmsg], SZ_LINE, matchfmt)
		    call pargstr (Memc[phottab])
		    call pargi (irow)

		    call synphotwarn (nomatch, Memc[errmsg])
		}
	    }

	    # Increment offset by the length of this table

	    offset = offset + nrow
	    call tbtclo (tp)
	}

	# Close files and release memory

	call clssyntab
	call freelist (optr)
	call freelist (sptr)
	call freelist (pptr)

	call sfree (sp)

end


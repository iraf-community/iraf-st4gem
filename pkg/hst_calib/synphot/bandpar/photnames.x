# PHOTNAMES -- Find which photometic names are included in the subset list

#* HISTORY *
#* B.Simon	06-Oct-94	original

procedure photnames (namelist, photlist, crosslist, maxphot, photvar, nphot)

char	namelist[ARB]	# i: complete list of photometric names
char	photlist[ARB]	# i: subset list of names
int	crosslist[ARB]	# i: list of names turned on when another turned on
int	maxphot		# i: number of names in complete list
int	photvar[ARB]	# o: set of flags indicating presence in sublist
int	nphot		# o: number of names in subset list
#--
bool	negate
int	ic, iphot, jphot
pointer	sp, word

string	nomatch  "Unknown photometric parameter"

bool	streq()
int	word_fetch(), word_match()

begin
	# Allocate memory for temporary string

	call smark (sp)
	call salloc (word, SZ_FNAME, TY_CHAR)

	negate = false
	call amovki (NO, photvar, maxphot)

	# Leading tilde is a negation sign, save info for later

	ic = 1
	while (photlist[ic] <= ' ')
	    ic = ic + 1

	if (photlist[ic] == '~') {
	    negate = true
	    ic = ic + 1
	}

	# Get each word in subset list, set appopriate flag

	iphot = 0
	while (word_fetch (photlist, ic, Memc[word], SZ_FNAME) > 0) {
	    iphot = iphot + 1
	    if (iphot > maxphot)
		break

	    call strlwr (Memc[word])
	    if (streq (Memc[word], "all")) {
		call amovki (YES, photvar, maxphot)
	    } else {
		jphot = word_match (Memc[word], namelist)
		if (jphot == 0)
		    call printerr_str (nomatch, Memc[word])

		photvar[jphot] = YES
	    }
	}

	# Set correlated flags according to their main flags

	do iphot = 1, maxphot {
	    jphot = crosslist[iphot]
	    if (jphot > 0) {
		if (photvar[jphot] == YES)
		    photvar[iphot] = YES
	    }
	}

	# Negate list of flags, if necessary

	if (negate) {
	    do iphot = 1, maxphot {
		if (photvar[iphot] == NO) {
		    photvar[iphot] = YES
		} else {
		    photvar[iphot] = NO
		}
	    }
	}

	# Count number of set flags

	nphot = 0
	do iphot = 1, maxphot {
	    if (photvar[iphot] == YES)
		nphot = nphot + 1
	}

	call sfree (sp)
end

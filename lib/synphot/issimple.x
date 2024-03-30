#* HISTORY *
#* B.Simon	20-Jun-94	do not treat numbers as simple bandpasses
#* B.Simon	16-Feb-94	rework logic of test
#* B.Simon	28-May-97	reworked logic again

# IS_SIMPLE -- Check if a bandpass command is simple

bool procedure is_simple (command)

char	command[ARB]	# i: command to be checked
#--
bool	simple
int	ic, nc, ntok, strict
pointer	sp, token, oldtok
real	rval

bool	isblank()
int	stridx(), tbtacc(), ctor(), syntok()

begin
	# Check command for a simple bandpass (no functions or filenames)

	call smark (sp)
	call salloc (token, SZ_FNAME, TY_CHAR)
	call salloc (oldtok, SZ_FNAME, TY_CHAR)

	# Bandpass is simple if ids are only separated by blanks or commas

	ic = 1
	ntok = 0
	strict = NO
	simple = true
	while (syntok (command, ic, strict, Memc[token], SZ_FNAME) > 0) {
	    ntok = ntok + 1
	    if (stridx (Memc[token], "()+-/*") > 0) {
		simple = false
		break
	    }

	    call strcpy (Memc[token], Memc[oldtok], SZ_FNAME)
	}

	# AND a single token is not a number

	if (ntok == 1 && simple) {
	    ic = 1
	    nc = ctor (Memc[oldtok], ic, rval)
	    simple = ! isblank (Memc[oldtok+ic-1])
	}

	# OR not a filename

	if (ntok == 1 && simple) {
	    simple = tbtacc (Memc[oldtok]) == NO
	}

	call sfree (sp)
	return (simple)
end

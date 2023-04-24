# CHKEXTEN -- Check a file names extension against a pattern
#
# B.Simon	16-Jan-91	Original

bool procedure chkexten (name, pattern)

char	name[ARB]	# i: file name
char	pattern[ARB]	# i: extension pattern
#--
bool	match
int	ic, junk
pointer	sp, cluster, ext

int	fnextn()

begin
	call smark (sp)
	call salloc (cluster, SZ_FNAME, TY_CHAR)
	call salloc (ext, SZ_FNAME, TY_CHAR)

	call imgcluster (name, Memc[cluster], SZ_FNAME)
	junk = fnextn (Memc[cluster], Memc[ext], SZ_FNAME)

	match = true
	for (ic = 1; match && pattern[ic] != EOS; ic = ic + 1) {
	    if (pattern[ic] != '?')
		match = pattern[ic] == Memc[ext+ic-1]
	}
	match = match && Memc[ext+ic-1] == EOS

	call sfree (sp)
	return (match)
end

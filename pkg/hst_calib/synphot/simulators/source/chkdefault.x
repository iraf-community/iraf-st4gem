#* HISTORY *
#* B.Simon	17-Feb-95	original

# CHKDEFAULT -- Replace the strings default and none with EOS

procedure chkdefault (str, maxch)

char	str[ARB]	# u: string to be checked
int	maxch		# i: maximum string length
#--
pointer	sp, temp
bool	streq()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (temp, maxch, TY_CHAR)

	# Convert string to lower case for simpler comparison

	call strcpy (str, Memc[temp], maxch)
	call strfix (Memc[temp])

	# If value is none or default, set to null string

	if (streq (Memc[temp], "default") || streq (Memc[temp], "none") ||
	    Memc[temp] == EOS) {
	    str[1] = EOS
	}

	call sfree (sp)

end

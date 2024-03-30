# MKMODESTR -- Make the observation mode string from an array of keywords

procedure mkmodestr (keyword, nkey, mode, maxch)

pointer	keyword[ARB]	# i: array of pointers to keyword strings
int	nkey		# i: number of keywords
char	mode[ARB]	# o: observation mode
int	maxch		# i: maximum length of observation mode
#--
int	ic, ikey
int	gstrcpy()

begin
	# The observation mode is a comma separated list of keywprds

	ic = 1
	do ikey = 1, nkey {
	    ic = ic + gstrcpy (Memc[keyword[ikey]], mode[ic], maxch-ic)
	    mode[ic] = ','
	    ic = ic + 1
	}

	if (ic > 1)
	    mode[ic-1] = EOS

end

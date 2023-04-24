# SORTFLAG -- Sort the list of flags
#
# B.Simon	11-Apr-91	First Code

procedure sortflag (flag, order, nflag)

int	flag[ARB]	# u: List of flag values
int	order[ARB]	# o: Order of flag in original list
int	nflag		# i: Number of flag values
#--
int	iflag, jflag
pointer	sp, temp

begin
	call smark (sp)
	call salloc (temp, nflag, TY_INT)

	call asrti (flag, Memi[temp], nflag)
	do iflag = 1, nflag {
	    for (jflag = 1; jflag <= nflag; jflag = jflag + 1) {
		if (Memi[temp+iflag-1] == flag[jflag])
		    break
	    }
	    order[iflag] = jflag
	}

	call amovi (Memi[temp], flag, nflag)
	call sfree (sp)
end

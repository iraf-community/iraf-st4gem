#* HISTORY *
#* B.Simon	21-Jul-94	original

# HPOVWRITE -- Read a table header parameter, possibly overwriting old value

procedure hpovwrite (tp, name, value, maxch)

pointer	tp		# i: table descriptor
char	name[ARB]	# i: header parameter name
char	value[ARB]	# u: header parameter value
int	maxch		# i: maximum length of value
#--
pointer	sp, oldval, newval, tmpval, errmsg

string	overwrite  "Task parameter overwritten by table header parameter"
string	msgformat  "%s = %s"

bool	strne()

begin
	call smark (sp)
	call salloc (oldval, maxch, TY_CHAR)
	call salloc (newval, maxch, TY_CHAR)
	call salloc (tmpval, maxch, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	ifnoerr (call tbhgtt (tp, name, Memc[tmpval], maxch)) {
	    call strcpy (Memc[tmpval], Memc[newval], maxch)
	    call strfix (Memc[newval])

	    if (Memc[newval] != EOS) {
		call strcpy (value, Memc[oldval], maxch)
		call strfix (Memc[oldval])

		if (strne (Memc[oldval], Memc[newval])) {
		    call strcpy (Memc[tmpval], value, maxch)

		    call sprintf (Memc[errmsg], SZ_LINE, msgformat)
		    call pargstr (name)
		    call pargstr (Memc[tmpval])

		    call synphotwarn (overwrite, Memc[errmsg])
		}
	    }
	}

	call sfree (sp)
end

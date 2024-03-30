# MODEWARN -- Conditionally print warning message about unused mode keywords

procedure modewarn (keystr)

char	keystr[ARB]			# i: list of unused keywords
int	value				# i: new value of flag
#--
int	warn				# print message about unused keywords?
data	warn  / YES /

string	notused  "Instrument mode keywords not used"

begin
	if (warn != NO)
	    call synphotwarn (notused, keystr)

	return

	# Set the variable controlling the mode warning message

	entry setmodewarn (value)

	warn = value
	return

end

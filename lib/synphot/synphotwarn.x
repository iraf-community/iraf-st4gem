# SYNPHOTWARN -- Write a message and take an error exit from synphot

procedure synphotwarn (errmsg, str)

char	errmsg[ARB]	# i: Error message
char	str[ARB]	# i: String to be appended to message
#--

begin
	call eprintf ("WARNING %s: %s\n")
	call pargstr (errmsg)
	call pargstr (str)

end


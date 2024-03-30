#* HISTORY *
#* B.Simon	28-Apr-94	original

# GETSYNVAR -- Get the value of a variable for a synphot expression

procedure getsynvar (ivar, value)

int	ivar		# i: variable number
real	value		# o: variable value
#--
real	variable[9]	# array of variable values
int	jvar

string	badivar  "Illegal variable number"
string	undefval "Value of variable is undefined"

begin
	if (ivar < 0 || ivar > 9)
	    call printerr_int (badivar, ivar)

	if (ivar == 0) {
	    call getvzero (value)
	} else {
	    value = variable[ivar]
	    if (IS_INDEFR(value))
		call printerr_int (undefval, ivar)
	}

	return

	# PUTSYNVAR -- Set the value of a variable

	entry putsynvar (ivar, value)

	if (ivar < 1 || ivar > 9)
	    call printerr_int (badivar, ivar)

	variable[ivar] = value
	return

	# UNDEFSYNVAR -- Set all variables to undefined

	entry undefsynvar ()

	call undefvzero
	do jvar = 1, 9
	    variable[jvar] = INDEFR
	return

end

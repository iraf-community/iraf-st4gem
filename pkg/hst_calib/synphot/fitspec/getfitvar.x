#* HISTORY *
#* B. Simon	05-Aug-94	original

# GETFIVAR -- Read least squares fit variables from parameter file

procedure getfitvar (var, nvar)

double	var[9]		# o: free variables in least squares fit
int	nvar		# o: number of variables
#--
int	ivar
pointer	sp, param

string	notcontig  "Fit variable is not contiguous with other variables"

double	clgetd()

begin
	call smark (sp)
	call salloc (param, SZ_FNAME, TY_CHAR)

	# Read the variable values

	for (ivar = 1; ivar <= 9; ivar = ivar + 1) {
	    call varname (ivar, Memc[param], SZ_FNAME)
	    var[ivar] = clgetd (Memc[param])
	}

	# Get number of variables that are not indef

	for (ivar = 1; ivar <= 9; ivar = ivar + 1) {
	    if (IS_INDEFD (var[ivar]))
		break
	}

	nvar = ivar - 1

	# Check for noncontiguous values

	for (ivar = nvar + 1; ivar <= 9; ivar = ivar + 1) {
	    if (! IS_INDEFD (var[ivar]))
		call printerr_int (notcontig, ivar)
	}

	call sfree (sp)

end

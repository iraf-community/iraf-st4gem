#* HISTORY *
#* B. Simon	05-Aug-94	original

# PUTFIVAR -- Write least squares fit variables to parameter file

procedure putfitvar (var, nvar)

double	var[9]		# i: free variables in least squares fit
int	nvar		# i: number of variables
#--
int	ivar
pointer	sp, param

begin
	call smark (sp)
	call salloc (param, SZ_FNAME, TY_CHAR)

	# Write the variable values. Break at the n-th variable.

	for (ivar = 1; ivar <= nvar; ivar = ivar + 1) {
	    call varname (ivar, Memc[param], SZ_FNAME)
	    call clputd (Memc[param], var[ivar])
	}

	call sfree (sp)
end

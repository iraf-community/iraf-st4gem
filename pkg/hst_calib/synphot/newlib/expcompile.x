#* HISTORY *
#* B.Simon	28-Apr-94	original
#* B.Simon	25-May-94	added call to exp_rewrite

# EXPCOMPILE -- Compile expression for syncalc

procedure expcompile (command, pcode, maxcode)

char	command[ARB]	# i: Command containing synphot expression
int	pcode[ARB]	# o: Pseudocode array for syncalc
int	maxcode		# o: Maximum length of pseudocode array
#--
pointer	sp, expr

begin
	call smark (sp)
	call salloc (expr, SZ_COMMAND, TY_CHAR)

	# Replace "naked" obsmodes with calls to the band() function

	call exp_rewrite (command, Memc[expr], SZ_COMMAND)
	call syncompile (Memc[expr], pcode, maxcode)

	call sfree (sp)
end


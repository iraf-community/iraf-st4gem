#* HISTORY *
#* B.Simon	07-Oct-94	original

# LENEXPR -- Get length of synphot expression

int procedure lenexpr (list)

pointer	list		# i: list descriptor for list of expressions
#--
int	len
pointer	sp, expr

int	strlen(), nxtlist()

begin
	# Allocate memory for temporary string

	call smark (sp)
	call salloc (expr, SZ_LINE, TY_CHAR)

	# Loop over each expression in list, saving the longest length

	len = 0
	while (nxtlist (list, Memc[expr], SZ_LINE) != EOF)
	    len = max (strlen (Memc[expr]), len)

	# Return longest length

	call sfree (sp)
	return (len)
end


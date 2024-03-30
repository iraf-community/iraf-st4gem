#* HISTORY *
#* B.Simon	13-May-94	Original
#* B.Simon	20-Jul-94	added qtlam

# FINDFUNC -- Find which function the user has requested

int procedure findfunc (form)

char	form[ARB]	# i: output form
#--
pointer	sp, lform
int	func

string	funclist "effstim,rmslam,fwhmlam,barlam,avglam,efflphot,efflerg"

int	word_match()

begin
	call smark (sp)
	call salloc (lform, SZ_FNAME, TY_CHAR)

	call strcpy (form, Memc[lform], SZ_FNAME)
	call strfix (Memc[lform])

	func = word_match (Memc[lform], funclist)	    

	call sfree (sp)
	return (func)

end

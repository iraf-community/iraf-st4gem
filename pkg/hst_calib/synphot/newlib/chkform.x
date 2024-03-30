#* HISTORY *
#* B.Simon	28-Apr-94	original

# CHKFORM -- Check the flux form for legality

procedure chkform (form)

char	form[ARB]	# i: output form
#--
pointer	sp, lform
int	word_match()

string	badform  "Unknown flux units"
string	formlist "photlam,counts,flam,fnu,photnu,jy,mjy,\
abmag,stmag,vegamag,obmag"		# list of valid flux units

begin
	call smark (sp)
	call salloc (lform, SZ_FNAME, TY_CHAR)

	call strcpy (form, Memc[lform], SZ_FNAME)
	call strfix (Memc[lform])

	if (word_match (Memc[lform], formlist) == 0)
	    call printerr_str (badform, form)

	call sfree (sp)
end

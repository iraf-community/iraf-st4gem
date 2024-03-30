include	"libsynphot.h"

#* HISTORY*
#* B.Simon	31-Oct-94	Remove error case for unrecognized units

int procedure is_count (units)

char	units[ARB]	# i: input flux units
#--
int	ftype, result, iscnt[11]
pointer	sp, form

data	iscnt  / NO, YES, 8 * NO, YES /

string	formlist FORMSTR

int	word_match()

begin
	call smark (sp)
	call salloc (form, SZ_FNAME, TY_CHAR)

	call strcpy (units, Memc[form], SZ_FNAME)
	call strfix (Memc[form])

	if (Memc[form] == EOS) {
	    result = NO

	} else {
	    ftype = word_match (Memc[form], formlist)
	    if (ftype == 0) {
		result = NO
	    } else {
		result = iscnt[ftype]
	    }
	}

	call sfree (sp)
	return (result)
end

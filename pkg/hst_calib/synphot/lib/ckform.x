define	SZ_FORM 10

# CKFORM -- Check that a given spectral or photometric form is valid

int procedure ckform( form )

char	form[ARB]	# i: Form to check
#--
int	strdic()

# Include the dictionary of forms
include	"dictionaries.h"

begin

	call strlwr( form )

	return( strdic( form, form, SZ_FORM, formdict ) )
end

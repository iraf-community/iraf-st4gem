#* HISTORY *
#* B.Simon	08-Apr-93	original
#* B.Simon	08-Jan-97	call is_image to do the work

include "imspec.h"

# GETFTYPE -- Get the type of the input file, even if there is no extension

int procedure getftype (file)

char	file[ARB]	# i: Name of file to be checked
#--
int	ftype
int	is_image(), tbtacc(), imaccess()

begin
	# Call is_image in the selector library to do the work

	switch (is_image (file)) {
	case ERR:
	    ftype = UNKNOWN

	case NO:
	    ftype = TABLE
	    if (tbtacc (file) == NO)
		ftype = NONE

	case YES:
	    ftype = IMAGE
	    if (imaccess (file, READ_ONLY) == NO)
		ftype = NONE
	}

	return (ftype)
end

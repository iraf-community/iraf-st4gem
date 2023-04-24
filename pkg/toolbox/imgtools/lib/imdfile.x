include <ctype.h>

# IMDFILE -- Create the name of the data file from the header file name
#
# B.Simon	02-Dec-88	First Code

procedure imdfile (hroot, droot, maxch)

char	hroot[ARB]	# i: Header file name
char	droot[ARB]	# o: Data file name
int	maxch		# i: Declared length of data file name
#--
int	junk
pointer	sp, extension, errmsg

string	badimgext  "Image extension not recognized (%s)"

bool	streq()
int	fnextn(), fnroot()

begin
	# Allocate storage space for strings

	call smark (sp)
	call salloc (extension, SZ_FNAME, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	# Extract extension from header name and check its legality

	if (fnextn (hroot, Memc[extension], SZ_FNAME) != 3) {
	    call sprintf (Memc[errmsg], SZ_LINE, badimgext)
		call pargstr (hroot)
	    call error (1, Memc[errmsg])
	}

	if (Memc[extension+2] != 'h') {
	    call sprintf (Memc[errmsg], SZ_LINE, badimgext)
		call pargstr (hroot)
	    call error (1, Memc[errmsg])
	}

	# Build data file name from header. Only two image
	# types are recognized -- OIF and STF

	junk = fnroot (hroot, droot, maxch)
	if (streq (Memc[extension], "imh")) {
	    call strcat (".pix", droot, maxch)
	} else {
	   Memc[extension+2] = 'd'
	   call strcat (".", droot, maxch)
	   call strcat (Memc[extension], droot, maxch)
	}

	call sfree (sp)
end

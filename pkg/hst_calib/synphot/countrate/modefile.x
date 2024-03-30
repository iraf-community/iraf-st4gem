#* HISTORY *
#* B.Simon	20-Mar-95	original
#* B.Simon	11-Jul-95	renamed, modified to call filemode

# MODEFILE -- Find a filename corresponding to to an observation mode 

bool procedure modefile (catalog, obsmode, output, maxch)

char	catalog[ARB]	# i: catalog name
char	obsmode[ARB]	# i: observation mode
char	output[ARB]	# o: output string
int	maxch		# i: length of output string
#--
int	jrow
pointer	sp, file, tp, cp

string	obscol   "OBSMODE"
string	filecol  "FILENAME"

pointer	tbtopn()

begin
	# Allocate memory for strings

	call smark (sp)
	call salloc (file, SZ_FNAME, TY_CHAR)

	# Open catalog

	call lastfile (catalog, Memc[file], SZ_FNAME)
	tp = tbtopn (Memc[file], READ_ONLY, 0)


	# Search for best match with observation mode

	call findmode (tp, obscol, obsmode, jrow)

	# Copy filename associated with obsmode to output string

	if (jrow == 0) {
	    output[1] = EOS

	} else {
	    call syncolptr (tp, filecol, 2, cp)
	    call tbegtt (tp, cp, jrow, output, maxch)
	}

	# Return flag which indicates success of match

	call sfree (sp)
	return (jrow != 0)

end

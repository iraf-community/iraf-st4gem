# LASTFILE -- Find the last file name which matches the template

# This procedure expands a file name template into a sorted list of file
# names and then selects the last name in the list. The file naming 
# conventions used in SDAS to generate unique names have the feature that 
# files generated later in time will also be later in alphabetical order.
# So this procedure can be used to retrieve the most recent version of a
# file of a given type.
#
# B.Simon	2-Oct-89	first code

procedure lastfile (template, fname, maxch)

char	template[ARB]		# i: File name template
char	fname[ARB]		# o: Last file which matches template
int	maxch			# i: Declared length of output string
#--
bool	found
int	index
pointer	list

int	fntrfnb(), fntlenb()
pointer	fntopnb()

begin
	found = false

	# Expand the file name template and retrieve the last name

	ifnoerr (list = fntopnb (template, YES)) {
	    index = fntlenb (list)
	    if (index > 0)
		found = fntrfnb (list, index, fname, maxch) > 0
	    call fntclsb (list)
	}

	# Copy the template to the file name if no match was found

	if (! found)
	    call strcpy (template, fname, maxch)

end

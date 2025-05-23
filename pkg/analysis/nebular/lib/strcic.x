# Copyright restrictions apply - see st4gem$copyright.stsdas

# Memory management
define	Dict_nc			Memc[dict_nc]
define	In_nc			Memc[in_nc]

#--------------------------------------------------------------------11 Apr 97--
.help strcic Mar96 nebular/lib
.ih
NAME
strcic -- Case-insensitive dictionary search
.ih
DESCRIPTION
Similar to strdic, except the string matching is case insensitive.
Both the dictionary and string to find are converted to lower case.
The string returned is how the string is found in the unmodified
dictionary.
.ih
ROUTINE DETAILS
.ls int procedure strcic (in_str, out_str, maxchars, dict)
.ls ARGUMENTS
.ls in_str[ARB] (input: char)
input string to find in dictionary
.le
.ls out_str[maxchars] (output: char)
string found in dictionary
.le
.ls maxchars (input: int)
maximum size of out_str
.le
.ls dict[ARB] (input: char)
dictionary
.le
.le
.le
.endhelp
#---------------------------------------------------------------------------
int procedure strcic (in_str, out_str, maxchars, dict)

char	in_str[ARB]		# I:  Input string to find in dictionary.
char	out_str[maxchars]	# O:  String found in dictionary.
int	maxchars		# I:  Maximum size of out_str.
char	dict[ARB]		# I:  Dictionary.

# Declarations
pointer	dict_nc			# Pointer to no-case dictionary.
int	dict_nc_sz		# Length of the no-case dictionary.
pointer	in_nc			# Pointer to no-case input string.
int	in_nc_sz		# Length of the no-case input string.
int	index			# Index of word in dictionary.
int	j			# Generic.
pointer	sp			# Stack pointer.
int	strlen()		# Length of string.
int	word_find()		# Retrive word from dictionary.
int	word_match()		# Find word in dictionary.

errchk	salloc, sfree, smark

begin
	call smark (sp)
	in_nc_sz = strlen (in_str)
	call salloc (in_nc, in_nc_sz, TY_CHAR)
	dict_nc_sz = strlen (dict)
	call salloc (dict_nc, dict_nc_sz, TY_CHAR)

	# Convert input string and dictionary to lower case.
	call strcpy (in_str, In_nc, in_nc_sz)
	call strlwr (In_nc)
	call strcpy (dict, Dict_nc, dict_nc_sz)
	call strlwr (Dict_nc)

	# Search the dictionary.
	index = word_match (In_nc, Dict_nc)

	# If something was found, fill the output string with the
	# correct form of the dictionary.
	if (index > 0)
	    j = word_find (index, dict, out_str, maxchars)

	# That's all folks.
	call sfree (sp)
	return (index)
end
#---------------------------------------------------------------------------
# End of strcic
#---------------------------------------------------------------------------

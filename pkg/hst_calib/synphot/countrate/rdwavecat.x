# RDWAVECAT -- Retrieve the name of a wavelength table from the catalog

bool procedure rdwavecat (wavecat, mode, wavetab, maxch)

char	wavecat[ARB]	# i: Wavelength catalog name
char	mode[ARB]	# i: Observation mode
char 	wavetab[ARB]	# o: Wavelength table name
int	maxch		# i: Max length of table name
#--
bool	match
pointer	sp, list

string	nowavetab  "Warning: cannot find wavelength table, using default\n"

bool	modefile()

begin
	# Allocate dynamic memory for temporary strings

	call smark (sp)
	call salloc (list, SZ_FNAME, TY_CHAR)

	# Extract list of obsmode keywords from mode string

	call getnaked (mode, Memc[list], SZ_FNAME)

	# Find mode string which matches input mode string

	match = modefile (wavecat, Memc[list], wavetab, SZ_FNAME)

	# Print warning message is wavelength file not found

	if (! match) {
	    call eprintf (nowavetab)
	    wavetab[1] = EOS
	}

	return (match)
end

include	"imspec.h"

#* HISTORY*
#* B.Simon	09-Apr-93	original

# RDTABFORM -- Read flux units from text table header

procedure rdtabform (table, inform, maxch)

char	table[ARB]	# i: table name
char	inform[ARB]	# o: spectral form string
int	maxch		# i: max length of spectral form string
#--
int	fd, ic
pointer	sp, line, word

string	funits   FLUXSTR

int	access(), open(), getline(), word_fetch(), word_match()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (line, SZ_LINE, TY_CHAR)
	call salloc (word, maxch, TY_CHAR)

	# If the table, is a text file search beginning commants for form

	inform[1] = EOS
	if (access (table, READ_ONLY, TEXT_FILE) == YES) {
	    fd = open (table, READ_ONLY, TEXT_FILE)

	    while (inform[1] == EOS && getline (fd, Memc[line]) != EOF) {
		if (Memc[line] != '#')
		    break

		# Check each word against legal forms for a match

		ic = 2
		while (word_fetch (Memc[line], ic, Memc[word], maxch) > 0) {
		    call strlwr (Memc[word])
		    if (word_match (Memc[word], funits) > 0) {
			call strcpy (Memc[word], inform, maxch)
			break
		    }
		}
	    }

	    call close (fd)
	}

	call sfree (sp)
end

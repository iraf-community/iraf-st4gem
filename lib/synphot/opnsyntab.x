include "libsynphot.h"

#* HISTORY *
#* B.Simon	23-Apr-93	original

# CLSSYNTAB -- Close cached tables

procedure clssyntab ()

#--
include "opnsyntab.com"

int	i

begin
	nextuse = 0
	initial = NO

	do i = 1, MAXCACHE {
	    if (tabptr[i] != NULL)
		call tbtclo (tabptr[i])

	    lastuse[i] = 0
	    tabptr[i] = NULL
	}

end

# ERRSYNTAB -- Error exit routine

procedure errsyntab (status)

int     status          # i: error status
#--

begin
	if (status != OK)
	    call clssyntab

end

# FREESYNTAB -- Remove a table from the cache

procedure freesyntab (table)

char	table[ARB]	# i: Table name
#--
include "opnsyntab.com"

int	i
pointer	sp, tabname

bool	streq()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (tabname, SZ_FNAME, TY_CHAR)

	# Find and close table

	do i = 1, MAXCACHE {
	    if (tabptr[i] == NULL)
		next

	    call tbtnam (tabptr[i], Memc[tabname], SZ_FNAME)
	    if (streq (table, Memc[tabname])) {
		call tbtclo (tabptr[i])
		tabptr[i] = NULL
	    }
	}

	call sfree (sp)
end

# INISYNTAB -- Initialize persistent variables used by opnsyntab

procedure inisyntab ()

#--
include "opnsyntab.com"

int	i
extern	errsyntab

begin
	if (initial == YES)
	    return

	nextuse = 0
	initial = YES
	do i = 1, MAXCACHE {
	    lastuse[i] = 0
	    tabptr[i] = NULL
	}

	call onerror (errsyntab)
end

# OPNSYNTAB -- Open synthetic photometry table or read from cache

pointer procedure opnsyntab (table)

char	table[ARB]	# i: Table name
#--
include "opnsyntab.com"

int	i, j
long	minuse
pointer	sp, tp, tabname

bool	streq()
pointer	tbtopn()

string	notinit "Synphot table cache not initialized"
string	nocache "Cannot open synphot table, no slot available in cache"
string	notopen "Cannot open synphot table, table open failed"

errchk	tbtopn, synphoterr

begin
	# Check initialization

	if (initial != YES)
	    call synphoterr (notinit, "opnsyntab")

	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (tabname, SZ_FNAME, TY_CHAR)

	# Search for table name in cache
	# At same time seach for empty slot in cache or
	# least recently accessed table

	nextuse = nextuse + 1
	minuse = nextuse
	j = 0

	do i = 1, MAXCACHE {
	    if (tabptr[i] == NULL) {
		minuse = 0
		j = i

	    } else {
		call syntabname (tabptr[i], Memc[tabname], SZ_FNAME)

		if (streq (table, Memc[tabname])) {
		    lastuse[i] = nextuse

		    call sfree (sp)
		    return (tabptr[i])
		}

		if (lastuse[i] < minuse) {
		    minuse = lastuse[i]
		    j = i
		}
	    }
	}

	# Error exit: could not find slot in cache

	if (j == 0)
	    call synphoterr (nocache, table)

	# Table not found in cache, open and add to cache

	if (tabptr[j] != NULL)
	    call tbtclo (tabptr[j])

	iferr {
	    tp = tbtopn (table, READ_ONLY, NULL)
	    lastuse[j] = nextuse
	    tabptr[j] = tp

	} then {
	    call synphoterr (notopen, table)
	    tp = NULL
	}

	call sfree (sp)
	return (tp)

end

include	"../grflist.h"

# OBSMODE -- Print the observation mode keywords from a graph table
#
# B.Simon	22-Jul-88	First Code

procedure obsmode ()

#--
pointer	path		# (Incomplete) observation mode path
pointer	grftable	# Graph table name

int	snode, ilist
pointer	sp, gp, keyname

string	nonode  "Observation mode path not found"
string	title1  "The complete list of observation mode keywords:\n"
string	title2  "The list of observation mode keywords for %s:\n"

bool	streq()
pointer	grf_open()

begin
	# Allocate dynamic memory for strings

	call smark (sp)
	call salloc (path, SZ_FNAME, TY_CHAR)
	call salloc (grftable, SZ_FNAME, TY_CHAR)

	# Read task parameters

	call clgnone ("path", Memc[path], SZ_FNAME)
	call clgstr ("grtbl", Memc[grftable], SZ_FNAME)

	# Read the instrument graph into a graph list structure

	gp = grf_open (Memc[grftable])
	call salloc (keyname, GRF_KEYLEN(gp), TY_CHAR)

	# Find starting component from observation mode

	if (Memc[path] == EOS) {
	    ilist = 0
	    call printf (title1)

	} else {
	    call startkey (gp, Memc[path], snode, Memc[keyname])

	    # Find starting node in graph table

	    for (ilist = 1; ilist <= GRF_SIZE(gp); ilist = ilist + 1) {
		if (GRF_INNODE(gp,ilist) == snode) {
		    if (streq (GRF_KEYNAME(gp,ilist), Memc[keyname])) {
			break
		    }
		}
	    }

	    if (ilist > GRF_SIZE(gp))
		call printerr_str (nonode, Memc[path])

	    # Print title for output:

	    call printf (title2)
	    call pargstr (Memc[path])
	}

	# Print keywords associated with descendants of start node

	call wrtmode (gp, ilist)

	# Free dynamic memory

	call grf_close (gp)
	call sfree (sp)

end

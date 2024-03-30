include "../grflist.h"

#* HISTORY *
#* B.Simon	03-Jan-95	Derived from getpath()
#* B.Simon	17-Mar-00	Reworked to handle multiple keywords better

# STARTPATH -- Find starting keyword from incomplete observation mode

procedure startpath (gp, complen, keylen, nmode, modelist, 
		     snode, keyname)

pointer	gp				# i: graph table descriptor
int	complen				# i: component name length
int	keylen				# i: keyword name length
int	nmode				# i: number of keywords in mode string
char	modelist[keylen,ARB]		# i: list of instrument mode keywords
int	snode				# o: starting node number
char	keyname[ARB]			# o: starting keyword name
#--
bool	uniq
int	irow, nrow, imode, umode
int	node, inode, onode, topnode, toplevel, newlevel
pointer	sp, used, keyword, compid, topkey, topid

string	notused   "Observation mode keyword not used"

bool	streq()

begin
	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (used, nmode, TY_INT)
	call salloc (keyword, keylen, TY_CHAR)
	call salloc (compid, complen, TY_CHAR)
	call salloc (topkey, keylen, TY_CHAR)
	call salloc (topid, complen, TY_CHAR)

	# Initialize variables used in searching graph

	umode = 0
	uniq = false
	node = GRF_INNODE(gp,1)
	call amovki (NO, Memi[used], nmode)
	call strcpy (GRF_COMPNAME(gp,1), Memc[compid], complen)

	topnode = 0
	toplevel = 0
	Memc[topid] = EOS

	# Initialize output variables

	call strcpy (GRF_KEYNAME(gp,1), keyname, keylen)
	snode = node

	nrow = GRF_SIZE(gp)
	for (irow = 1; irow <= nrow; irow = irow + 1) {

	    # Find the next row in the graph whose input node is
	    # equal to the output node of the previous row in the graph

	    inode = GRF_INNODE(gp,irow)

	    if (inode == node) {
		# Read the keyword, component name and 
		# output node from the table

		call strcpy (GRF_KEYNAME(gp,irow), Memc[keyword], keylen)
		call strcpy (GRF_COMPNAME(gp,irow), Memc[compid], complen)

		onode = GRF_OUTNODE(gp,irow)

		# Find the best continuation of the current path
		# level = 3 means that we have matched a keyword
		# level = 2 means that we have found the default path
		# level = 1 means that we have found a possible continuation

		newlevel = 1

		# Only traverse default nodes if there are unused keywords

		if (streq (Memc[keyword], "default") && umode < nmode)
		    newlevel = 2

		do imode = 1, nmode {
		    if (streq (Memc[keyword], modelist[1,imode])) {
			if (Memi[used+imode-1] == NO) {
			    Memi[used+imode-1] = YES

			    # The case of only a single keyword is handled
			    # as a special case because of the way the graph
			    # table is structured. Otherwise, every node would
			    # be returned when only the instrument is specified

			    if (nmode > 1)
				umode = umode + 1
			}

			call strcpy (Memc[keyword], keyname, keylen)
			snode = node
			newlevel = 3
		    }
		}

		# Save information if this is the best match so far
		# If more than one row has the same level, uniq is false

		if (newlevel > toplevel) {
		    call strcpy (Memc[keyword], Memc[topkey], keylen)
		    call strcpy (Memc[compid], Memc[topid], complen)
		    toplevel = newlevel
		    topnode = onode
		    uniq = true

		} else if (newlevel == toplevel) {
		    uniq = streq (Memc[topid], Memc[compid]) &&
			   topnode == onode
		}
	    }

	    if ((inode != node || irow == nrow) && toplevel > 0) {
		# Exit loop if no unique continuation can be found

		if (! uniq)
		    break

		Memc[topid] = EOS
		node = topnode
		toplevel = 0
		topnode = 0

		irow = irow - 1
	    }
	}

	# Send error message if not all keywords used

	if (umode < nmode) {
	    do imode = 1, nmode {
		if (Memi[used+imode-1] == NO)
		    call printerr_str (notused, modelist[1,imode])
	    }
	}

	call sfree (sp)
end

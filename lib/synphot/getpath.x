include "libsynphot.h"

# GETPATH -- Trace path of optical components through the graph table

procedure getpath (verbose, grf, grf_index, mxlist, nmode, ncomp, 
		   modelist, modenum, complist)

bool	verbose				# i: diagnostic message switch
pointer	grf				# i: graph table structure
int	grf_index			# i: selects the column in the table to return
int	mxlist				# i: maximum number of keywords in list
int	nmode				# i: number of keywords in mode string
int	ncomp				# o: number of instrument components
char	modelist[SZ_KEYWRD,ARB]		# i: list of instrument mode keywords
int	modenum[ARB]			# o: component number of mode keyword
char	complist[SZ_COMPID,ARB]		# o: list of instrument components
#--
bool	uniq
int	ic, irow, jrow, nrow, imode
int	node, inode, onode, topnode, toplevel, newlevel
pointer	sp, keylist, keyword, compid, topid

string	header    "\nList of component names:\n"
string	bodyfmt   "\t%s\n"
string	notuniq   "Mode incomplete. Continuations"
string	overcomp  "List of components overflows array. Last component"
string	backnode  "Output node in graph table points backward. Component name"

bool	streq()
int	gstrcpy(), sizegraf()

errchk	synphoterr, rdgrafi, rdgraft

begin


	# Allocate memory for temporary strings

	call smark (sp)
	call salloc (keylist, SZ_COMMAND, TY_CHAR)
	call salloc (keyword, SZ_KEYWRD, TY_CHAR)
	call salloc (compid, SZ_COMPID, TY_CHAR)
	call salloc (topid, SZ_COMPID, TY_CHAR)

	# Initialize output variables

	ncomp = 0
	call amovki (0, modenum, nmode)

	# Print header for list of components

	if (verbose)
	    call eprintf (header)

	# Initialize variables used in searching graph

	jrow = 0
	uniq = false
	call rdgrafi (grf, GRF_INNODE, 1, node)

	topnode = 0
	toplevel = 0
	Memc[topid] = EOS

	nrow = sizegraf (grf)
	for (irow = 1; irow <= nrow; irow = irow + 1) {

	    # Find the next row in the graph whose input node is
	    # equal to the output node of the previous row in the graph

	    call rdgrafi (grf, GRF_INNODE, irow, inode)

	    if (inode == node) {
		# Read the keyword, component name and 
		# output node from the table

		call rdgraft (grf, GRF_KEYWRD, irow, Memc[keyword], SZ_KEYWRD)
		call strfix (Memc[keyword])

		call rdgraft (grf, grf_index, irow, Memc[compid], SZ_COMPID)
		call strfix (Memc[compid])

		call rdgrafi (grf, GRF_OUTNODE, irow, onode)

		# Find the best continuation of the current path
		# level = 3 means that we have matched a keyword
		# level = 2 means that we have found the default path
		# level = 1 means that we have found a possible continuation

		newlevel = 1

		if (streq (Memc[keyword], "default"))
		    newlevel = 2

		do imode = 1, nmode {
		    if (streq (Memc[keyword], modelist[1,imode])) {
			modenum[imode] = ncomp + 1
			newlevel = 3
			break
		    }
		}

		# Save information if this is the best match so far
		# If more than one row has the same level, uniq is false

		if (newlevel > toplevel) {
		    call strcpy (Memc[compid], Memc[topid], SZ_COMPID)
		    toplevel = newlevel
		    topnode = onode
		    uniq = true

		} else if (newlevel == toplevel) {
		    uniq = streq (Memc[topid], Memc[compid]) &&
			   topnode == onode
		}

		# Jrow tracks the first row with the proper input node

		if (jrow == 0)
		    jrow = irow
	    }

	    if ((inode != node || irow == nrow) && toplevel > 0) {
		# Print the list of possible continuations
		# if no unique continuation can be found

		if (! uniq) {
		    ic = 0
		    while (jrow < irow) {
			call rdgraft (grf, GRF_KEYWRD, jrow, 
				      Memc[keyword], SZ_KEYWRD)
			call strfix (Memc[keyword])

			ic = ic + gstrcpy (Memc[keyword], Memc[keylist+ic],
					   SZ_LINE-(ic+1))

			Memc[keylist+ic] = ' '
			ic = ic + 1
			jrow = jrow + 1
		    }

		    Memc[keylist+ic] = EOS
		    call synphoterr (notuniq, Memc[keylist])
		}

		# Skip clear components, add other components
		# to the list of component names

		if (streq (Memc[topid], "clear")) {
		    do imode = 1, nmode {
			if (modenum[imode] == ncomp + 1)
			    modenum[imode] = -1
		    }

		} else {
		    ncomp = ncomp + 1
		    if (ncomp > mxlist)
			call synphoterr (overcomp, Memc[topid])
		    call strcpy (Memc[topid], complist[1,ncomp], SZ_COMPID)

		    if (verbose) {
			call eprintf (bodyfmt)
			call pargstr (Memc[topid])
		    }
		}

		if (topnode <= node) {
		    call printf( "topnode, node: %d %d \n")
		    call pargi(topnode)
		    call pargi(node)
		    call flush(STDOUT)
		    call synphoterr (backnode, Memc[topid])
		}
		node = topnode

		Memc[topid] = EOS
		toplevel = 0
		topnode = 0

		jrow = 0
		irow = irow - 1
	    }
	}

	# If the list is empty, put a clear component on this list

	if (ncomp == 0) {
	    ncomp = 1
	    call strcpy ("clear", complist[1,ncomp], SZ_COMPID)
	}

	call sfree (sp)
end

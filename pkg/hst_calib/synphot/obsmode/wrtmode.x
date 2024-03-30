include "../grflist.h"

#* HISTORY *
#* B.Simon	04-Jan-95	original

# WRTMODE -- Write the keywords associated with an observation mode

procedure wrtmode (gp, root)

pointer	gp		# i: graph table descriptor
int	root		# i: index to start node
#--
int	nlist, bufsize, top, nkey, index, ilist, node, jlist, klist
pointer	sp, inode, keyidx, keybuf, hash

bool	strne(), add_hash()
int	gstrcpy()
pointer	opn_hash()

begin
	# Allocate temporary arrays

	nlist = GRF_SIZE(gp)
	bufsize = nlist * (GRF_KEYLEN(gp) + 1)

	call smark (sp)
	call salloc (inode, nlist, TY_INT)
	call salloc (keyidx, nlist, TY_INT)
	call salloc (keybuf, bufsize, TY_CHAR)

	# Initialize variables that track keyword strings

	nkey = 0
	index = 1
	hash = opn_hash (nlist, GRF_KEYLEN(gp))

	# Initialize list of nodes that are descended from the root node

	top = 0
	if (root == 0) {
	    Memi[inode] = GRF_INNODE(gp,1)
	} else {
	    Memi[inode] = GRF_OUTNODE(gp,root)
	}

	# Print each non-default keyword that is descended from the root

	for (ilist = root + 1; ilist <= nlist; ilist = ilist + 1) {
	    # Exit if no pending nodes in list

	    if (top < 0)
		break

	    node = GRF_INNODE(gp,ilist)
	    if (node == Memi[inode+top]) {
		# Add keyword to list if not default

		if (strne(GRF_KEYNAME(gp,ilist), "default")) {
		    if (add_hash (hash, GRF_KEYNAME(gp,ilist))) {
			Memi[keyidx+nkey] = index
			nkey = nkey + 1

			index = gstrcpy (GRF_KEYNAME(gp,ilist), 
					 Memc[keybuf+index-1], 
					 bufsize - index) +
					 index + 1
		    }
		}

		# Add output node to list of descendants

		for (jlist = 0; jlist <= top; jlist = jlist + 1) {
		    if (Memi[inode+jlist] <= GRF_OUTNODE(gp,ilist))
			break
		}

		# Insert in proper sort order in list

		if (jlist <= top) {
		    if (Memi[inode+jlist] != GRF_OUTNODE(gp,ilist)) {
			do klist = top, jlist, -1
			    Memi[inode+klist+1] = Memi[inode+klist]

			Memi[inode+jlist] = GRF_OUTNODE(gp,ilist)
			top = top + 1
		    }
		}

	    } else if (node > Memi[inode+top]) {
		# Write list of keywords to STDOUT

		call wrtkeys (Memi[keyidx], Memc[keybuf], nkey, GRF_KEYLEN(gp))

		# Reinitialize list of keywords and drop node from list

		nkey = 0
		index = 1
		top = top - 1
		ilist = ilist - 1
	    }
	}

	call wrtkeys (Memi[keyidx], Memc[keybuf], nkey, GRF_KEYLEN(gp))

	call free_hash(hash)
	call sfree (sp)
end
